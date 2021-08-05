#' @export
SBC_results <- function(stats, fits, errors = rep(list(NULL), length(fits))) {
  validate_SBC_results(
    structure(list(stats = stats, fits = fits, errors = errors), class = "SBC_results")
  )
}

#' @export
validate_SBC_results <- function(x) {
  stopifnot(is.list(x))
  stopifnot(inherits(x, "SBC_results"))
  if(!is.data.frame(x$stats)) {
    stop("SBC_datasets object has to have a 'stats' field of type data.frame")
  }

  if(!is.list(x$fits)) {
    stop("SBC_datasets object has to have a 'fits' field of type list")
  }

  if(!is.list(x$errors)) {
    stop("SBC_datasets object has to have an 'errors' field of type list")
  }

  if(!is.numeric(x$stats$run_id) || min(x$stats$run_id) <= 0) {
    stop("The run_id column of stats needs to be number > 0")
  }

  if(length(unique(x$stats$run_id)) != length(x$fits)) {
    stop("Needs equal no. of stats and fits")
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
  errors_list <- purrr::map(args, function(x) x$errors)

  # Ensure unique run_ids
  max_ids <- as.numeric(purrr::map(stats_list, function(x) max(x$run_id)))
  shifts <- c(0, max_ids[1:(length(max_ids)) - 1])

  stats_list <- purrr::map2(stats_list, shifts, function(x, shift) dplyr::mutate(x, run_id = run_id + shift))

  SBC_results(do.call(rbind, stats_list),
                   do.call(c, fits_list),
              do.call(c, errors_list))
}

#' @export
compute_results <- function(datasets, backend,
                            cores_per_fit = default_cores_per_fit(length(datasets)),
                            keep_fits = TRUE,
                            thin_ranks = 1,
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
  n_errors <- 0
  max_errors_to_show <- 5
  for(i in 1:length(datasets)) {
    if(is.null(results_raw[[i]]$error)) {
      fits[[i]] <- results_raw[[i]]$fit
      stats_list[[i]] <- results_raw[[i]]$stats
      stats_list[[i]]$run_id <- i
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

  if(length(unique(stats$max_rank)) != 1) {
    warning("Differening max_rank across fits")
  }
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

  SBC_results(stats = stats, fits = fits, errors = errors)
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
    if(!keep_fit) {
      fit <- NULL
    }
    c(list(fit = fit, stats = stats, error = NULL))
  }, error = function(e) { list(fit = NULL, stats = NULL, error = e) })
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

