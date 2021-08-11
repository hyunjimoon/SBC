#' @export
SBC_results <- function(stats,
                        fits,
                        fit_diagnostics = NULL,
                        errors = rep(list(NULL), length(fits))) {
  param_diagnostics <-  tryCatch(
    { dplyr::summarise(dplyr::group_by(stats, run_id),
                                n_params = dplyr::n(),
                                max_rhat = max(rhat),
                                min_ess = min(ess_bulk),
                                min_ess_to_rank = min(ess_bulk / max_rank)) },
    error = function(e) { warning("Error when computing param diagnostics. ", e)})
  validate_SBC_results(
    structure(list(stats = stats, fits = fits, fit_diagnostics = fit_diagnostics,
                   param_diagnostics = param_diagnostics, errors = errors), class = "SBC_results")
  )
}

#' @export
validate_SBC_results <- function(x) {
  stopifnot(is.list(x))
  stopifnot(inherits(x, "SBC_results"))
  if(!is.data.frame(x$stats)) {
    stop("SBC_results object has to have a 'stats' field of type data.frame")
  }

  if(!is.list(x$fits)) {
    stop("SBC_results object has to have a 'fits' field of type list")
  }

  if(!is.null(x$fit_diagnostics) && !is.data.frame(x$fit_diagnostics)) {
    stop("If the SBC_results object has a 'fit_diagnostics' field, it has to inherit from data.frame")
  }

  if(!is.data.frame(x$param_diagnostics)) {
    stop("If the SBC_results object has a 'param_diagnostics' field, it has to inherit from data.frame")
  }


  if(!is.list(x$errors)) {
    stop("SBC_results object has to have an 'errors' field of type list")
  }

  if(nrow(x$stats) > 0) {
    if(!is.numeric(x$stats$run_id)) {
      stop("The run_id column of stats needs to be a number.")
    }


    if(min(x$stats$run_id) < 1 || max(x$stats$run_id > length(x$fits))) {
      stop("stats$run_id values must be between 1 and number of fits")
    }
  }

  if(!is.null(x$fit_diagnostics) && nrow(x$fit_diagnostics) > 0) {
    if(!is.numeric(x$fit_diagnostics$run_id)) {
      stop("The run_id column of 'fit_diagnostics' needs to be a number.")
    }


    if(min(x$fit_diagnostics$run_id) < 1 || max(x$fit_diagnostics$run_id > length(x$fits))) {
      stop("fit_diagnostics$run_id values must be between 1 and number of fits")
    }
  }

  if(nrow(x$param_diagnostics) > 0) {
    if(!is.numeric(x$param_diagnostics$run_id)) {
      stop("The run_id column of 'param_diagnostics' needs to be a number.")
    }


    if(min(x$param_diagnostics$run_id) < 1 || max(x$param_diagnostics$run_id > length(x$fits))) {
      stop("param_diagnostics$run_id values must be between 1 and number of fits")
    }
  }


  if(length(x$fits) != length(x$errors)) {
    stop("Needs equal no. of fits and errors")
  }

  #TODO check identical par names
  x
}


#' @export
bind_results <- function(...) {
  args <- list(...)

  purrr::walk(args, validate_SBC_results)


  stats_list <- purrr::map(args, function(x) x$stats)
  fits_list <- purrr::map(args, function(x) x$fits)
  fit_diagnostics_list <- purrr::map(args, function(x) x$fit_diagnostics)
  errors_list <- purrr::map(args, function(x) x$errors)

  # Ensure unique run_ids
  max_ids <- as.numeric(purrr::map(stats_list, function(x) max(x$run_id)))
  shifts <- c(0, max_ids[1:(length(max_ids)) - 1])

  shift_run_id <- function(x, shift) {
    if(is.null(x)) {
      x
    } else {
      dplyr::mutate(x, run_id = run_id + shift)
    }
  }

  stats_list <- purrr::map2(stats_list, shifts, shift_run_id)
  fit_diagnostics_list <- purrr::map2(fit_diagnostics_list, shifts, shift_run_id)

  SBC_results(stats = do.call(rbind, stats_list),
              fits = do.call(c, fits_list),
              fit_diagnostics = do.call(rbind, fit_diagnostics_list),
              errors =  do.call(c, errors_list))
}

#' @export
compute_results <- function(datasets, backend,
                            cores_per_fit = default_cores_per_fit(length(datasets)),
                            keep_fits = TRUE,
                            thin_ranks = 10,
                            chunk_size = default_chunk_size(length(datasets)) ) {
  stopifnot(length(datasets) > 0)

  # Create combined data for computation
  params_and_generated_list <- list()
  for(i in 1:length(datasets)) {
    params_and_generated_list[[i]] <- list(
      parameters = posterior::subset_draws(datasets$parameters,
                                      draw = i),
      generated = datasets$generated[[i]]
    )
  }

  results_raw <- future.apply::future_lapply(
    params_and_generated_list, SBC:::compute_results_single,
    backend = backend, cores = cores_per_fit,
    keep_fit = keep_fits, thin_ranks = thin_ranks,
    future.seed = TRUE,
    future.chunk.size = chunk_size)


  # Combine, check and summarise
  fits <- rep(list(NULL), length(datasets))
  errors <- rep(list(NULL), length(datasets))
  stats_list <- list()
  fit_diagnostics_list <- list()
  n_errors <- 0
  max_errors_to_show <- 5
  for(i in 1:length(datasets)) {
    if(is.null(results_raw[[i]]$error)) {
      fits[[i]] <- results_raw[[i]]$fit
      stats_list[[i]] <- results_raw[[i]]$stats
      stats_list[[i]]$run_id <- i
      fit_diagnostics_list[[i]] <- results_raw[[i]]$fit_diagnostics
      fit_diagnostics_list[[i]]$run_id <- i
    }
    else {
      if(n_errors < max_errors_to_show) {
        warning("Dataset ", i, " resulted in error when fitting.\n")
        message(results_raw[[i]]$error, "\n")
      } else if(n_errors == max_errors_to_show) {
        warning("Too many datasets produced errors. Further error messages not shown.\n")
      }
      n_errors <- n_errors + 1
      errors[[i]] <- results_raw[[i]]$error
    }
  }

  if(n_errors == length(datasets)) {
    stop("All datasets produced error when fitting")
  } else if(n_errors > 0) {
    warning("Total of ", n_errors, " datasets produced errors.")
  }

  stats <- do.call(rbind, stats_list)

  unique_max_ranks <- unique(stats$max_rank)
  if(length(unique_max_ranks) != 1) {
    warning("Differening max_rank across fits")
  }

  if(min(unique_max_ranks) < 50) {
    warning("Ranks were computed from fewer than 50 samples, the SBC checks will have low ",
            "precision.\nYou may need to increase the number of samples from the backend and make sure that ",
            "the combination of thinning in the backend and `thin_ranks` is sensible.\n",
            "Currently thin_ranks = ", thin_ranks, ".")

  }

  fit_diagnostics <- do.call(rbind, fit_diagnostics_list)

  all_vars <- dplyr::summarise(
    dplyr::group_by(stats, run_id),
    all_vars = paste0(variable, collapse = ","), .groups = "drop")
  if(length(unique(all_vars$all_vars)) != 1) {
    warning("Not all fits share the same variables")
  }

  missing_vars <- setdiff(posterior::variables(datasets$parameters), stats$variable)
  if(length(missing_vars) > 0) {
    warning("Some variables missing in fits: ", paste0(missing_vars, collapse = ", "))

  }

  SBC_results(stats = stats, fits = fits, fit_diagnostics = fit_diagnostics, errors = errors)
}

#' @export
default_cores_per_fit <- function(n_fits, total_cores = future::availableCores(),
                                  chunk_size = default_chunk_size(n_fits)) {
  if(inherits(future::plan(), "sequential")) {
    total_cores
  } else if(2 * (n_fits / chunk_size) <= total_cores) {
    floor(total_cores / (n_fits / chunk_size))
  } else {
    1
  }
}

#' @export
default_chunk_size <- function(n_fits, n_workers = future::nbrOfWorkers()) {
  guess <- if(is.infinite(n_workers)) {
    1
  } else {
    n_fits / n_workers
  }
  max(guess, getOption("SBC.min_chunk_size", 1))
}



compute_results_single <- function(params_and_generated, backend, cores,
                                   keep_fit, thin_ranks) {
  tryCatch({
    parameters <- params_and_generated$parameters
    generated <- params_and_generated$generated
    fit <- SBC_fit(backend, generated, cores = cores)
    stats <- statistics_from_single_fit(fit, parameters = parameters, thin_ranks = thin_ranks)
    fit_diagnostics <- SBC_fit_to_diagnostics(fit)
    if(!keep_fit) {
      fit <- NULL
    }
    c(list(fit = fit, stats = stats, fit_diagnostics = fit_diagnostics, error = NULL))
  }, error = function(e) { list(fit = NULL, stats = NULL, fit_diagnostics = NULL, error = e) })
}

#' @export
statistics_from_single_fit <- function(fit, parameters, thin_ranks) {
  fit_matrix <- SBC_fit_to_draws_matrix(fit)

  shared_vars <- intersect(posterior::variables(parameters),
                           posterior::variables(fit_matrix))


  # Make sure the order of variables matches
  parameters <- posterior::subset_draws(parameters, variable = shared_vars)


  fit_matrix <- posterior::subset_draws(fit_matrix, variable = shared_vars)

  fit_thinned <- posterior::thin_draws(fit_matrix, thin_ranks)


  stats <- posterior::summarise_draws(fit_matrix)
  stats$simulated_value <- as.numeric(parameters)

  ranks <- calculate_ranks_draws_matrix(parameters, fit_thinned)
  if(!identical(stats$variable, names(ranks))) {
    stop("A naming conflict")
  }
  stats$rank <- ranks
  stats$max_rank <- attr(ranks, "max_rank")
  stats$z_score <- (stats$simulated_value - stats$mean) / stats$sd

  stats <- dplyr::select(
    stats, simulated_value, rank, z_score, tidyselect::everything())

  stats
}

#' Discrete uniform distribution allowing for varying lower and upper bounds.
#'
#' Internal, should not be exported.
#' Based on https://stats.stackexchange.com/a/3939/73129
rdunif <- function(n, a, b) {
  ceiling(runif(n, min = a - 1, max= b))
}

#' @param dm draws_matrix of the fit (assumed to be already thinned if that was necessary)
#' @export
calculate_ranks_draws_matrix <- function(params, dm) {
  #TODO validate input
  max_rank <- posterior::ndraws(dm)

  less_matrix <- sweep(dm, MARGIN = 2, STATS = params, FUN = "<")
  rank_min <- colSums(less_matrix)

  # When there are ties (e.g. for discrete parameters), the rank is currently drawn stochastically
  # among the ties
  equal_matrix <- sweep(dm, MARGIN = 2, STATS = params, FUN = "==")
  rank_range <- colSums(equal_matrix)

  ranks <- rank_min + rdunif(posterior::nvariables(dm), a = 0, b = rank_range)

  attr(ranks, "max_rank") <- max_rank
  ranks
}

#' @export
calculate_sds_draws_matrix <- function(dm) {
  #TODO: validate input
  apply(dm, MARGIN = 2, FUN = sd)
}

#' @export
summary.SBC_results <- function(x) {
  summ <- list(
    n_fits = length(x$fits),
    n_errors = sum(!purrr::map_lgl(x$errors, is.null)),
    n_high_rhat = sum(x$param_diagnostics$max_rhat > 1.01),
    n_low_ess_to_rank = sum(x$param_diagnostics$min_ess_to_rank < 0.5))
  if(!is.null(x$fit_diagnostics)) {
    summ$fit_diagnostics <- summary(x$fit_diagnostics)
  } else {
    summ$fit_diagnostics <- NULL
  }
  structure(
    summ,
    class = "SBC_results_summary"
  )
}

#' @export
print.SBC_results_summary <- function(x) {
  cat("SBC_results with ", x$n_fits, " total fits.\n")
  if(x$n_errors > 0) {
    cat(" - ", x$n_errors, " (", round(100 * x$n_errors / x$n_fits), "%) fits resulted in an error.\n", sep = "")
  } else {
    cat(" - No fits had errors.\n")
  }
  if(x$n_high_rhat > 0) {
    cat(" - ", x$n_high_rhat, " (", round(100 * x$n_high_rhat / x$n_fits), "%) fits had at least one Rhat > 1.01.\n", sep = "")
  } else {
    cat(" - No fits had Rhat > 1.01.\n")
  }
  if(x$n_low_ess_to_rank > 0) {
    cat(" - ", x$n_low_ess_to_rank, " (", round(100 * x$n_low_ess_to_rank / x$n_fits), "%) fits had bulk ESS < ",
    "half of the maximum rank, potentially skewing the rank statistics.\n   Consider increasing `thin_ranks` or number of posterior samples and recomputing.", sep = "")
  } else {
    cat(" - All fits had bulk ESS > half of the maximum rank.\n")
  }

  if(!is.null(x$fit_diagnostics)) {
    print(x$fit_diagnostics)
  }

  invisible(x)
}
