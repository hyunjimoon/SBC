#' SBC_results objects.
#'
#'
#' The `SBC_results` contains the following fields:
#'   - `$stats` statistics for all parameters and fits (one row per parameter-fit combination)
#'   - `$fits`  the raw fits (unless `keep_fits = FALSE`) or `NULL` if the fit failed
#'   - `$errors` error messages that caused fit failures
#'   - `$outputs`, `$messages`, `$warnings` the outputs/messages/warnings written by fits
#'   - `$default_diagnostics` a data frame of default convergence/correctness diagnostics (one row per fit)
#'   - `$backend_diagnostics` a data frame of backend-specific diagnostics (one row per fit)
#' @export
SBC_results <- function(stats,
                        fits,
                        backend_diagnostics,
                        default_diagnostics,
                        outputs,
                        messages,
                        warnings,
                        errors) {
  validate_SBC_results(
    structure(list(stats = stats, fits = fits, backend_diagnostics = backend_diagnostics,
                   outputs = outputs, messages = messages, warnings = warnings,
                   default_diagnostics = default_diagnostics, errors = errors), class = "SBC_results")
  )
}

compute_default_diagnostics <- function(stats) {
  dplyr::summarise(dplyr::group_by(stats, dataset_id),
                   n_params = dplyr::n(),
                   max_rhat = max(c(-Inf, rhat)),
                   min_ess_bulk = min(c(Inf, ess_bulk)),
                   min_ess_tail = min(c(Inf, ess_tail)),
                   min_ess_to_rank = min(c(Inf, ess_tail / max_rank)))
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

  if(!is.null(x$backend_diagnostics) && !is.data.frame(x$backend_diagnostics)) {
    stop("If the SBC_results object has a 'backend_diagnostics' field, it has to inherit from data.frame")
  }

  if(!is.data.frame(x$default_diagnostics)) {
    stop("If the SBC_results object has a 'default_diagnostics' field, it has to inherit from data.frame")
  }


  if(!is.list(x$errors)) {
    stop("SBC_results object has to have an 'errors' field of type list")
  }

  if(nrow(x$stats) > 0) {
    if(!is.numeric(x$stats$dataset_id)) {
      stop("The dataset_id column of stats needs to be a number.")
    }


    if(min(x$stats$dataset_id) < 1 || max(x$stats$dataset_id) > length(x$fits)) {
      stop("stats$dataset_id values must be between 1 and number of fits")
    }
  }

  if(!is.null(x$outputs)) {
    if(!is.list(x$outputs) || length(x$outputs) != length(x$fits)) {
      stop("outputs can only be a list of the same length as fits")
    }
  }

  if(!is.null(x$messages)) {
    if(!is.list(x$messages) || length(x$messages) != length(x$fits)) {
      stop("messages can only be a list of the same length as fits")
    }
  }

  if(!is.null(x$warnings)) {
    if(!is.list(x$warnings) || length(x$warnings) != length(x$fits)) {
      stop("warnings can only be a list of the same length as fits")
    }
  }

  if(!is.null(x$backend_diagnostics) && nrow(x$backend_diagnostics) > 0) {
    if(!is.numeric(x$backend_diagnostics$dataset_id)) {
      stop("The dataset_id column of 'backend_diagnostics' needs to be a number.")
    }


    if(min(x$backend_diagnostics$dataset_id) < 1 || max(x$backend_diagnostics$dataset_id > length(x$fits))) {
      stop("backend_diagnostics$dataset_id values must be between 1 and number of fits")
    }
  }

  if(nrow(x$default_diagnostics) > 0) {
    if(!is.numeric(x$default_diagnostics$dataset_id)) {
      stop("The dataset_id column of 'default_diagnostics' needs to be a number.")
    }


    if(min(x$default_diagnostics$dataset_id) < 1 || max(x$default_diagnostics$dataset_id > length(x$fits))) {
      stop("default_diagnostics$dataset_id values must be between 1 and number of fits")
    }
  }


  if(length(x$fits) != length(x$errors)) {
    stop("Needs equal no. of fits and errors")
  }

  #TODO check identical par names
  x
}

#' Combine multiple SBC results together.
#'
#' Primarily useful for iteratively adding more datasets to your SBC check.
#'
#' An example usage can be found in the `small_model_workflow` vignette.
#' @param ... objects of type `SBC_results` to be combined.
#' @export
bind_results <- function(...) {
  args <- list(...)

  purrr::walk(args, validate_SBC_results)


  stats_list <- purrr::map(args, function(x) x$stats)
  fits_list <- purrr::map(args, function(x) x$fits)
  backend_diagnostics_list <- purrr::map(args, function(x) x$backend_diagnostics)
  default_diagnostics_list <- purrr::map(args, function(x) x$default_diagnostics)
  errors_list <- purrr::map(args, function(x) x$errors)
  messages_list <- purrr::map(args, function(x) x$messages)
  warnings_list <- purrr::map(args, function(x) x$warnings)
  outputs_list <- purrr::map(args, function(x) x$outputs)

  # Ensure unique dataset_ids
  max_ids <- as.numeric(purrr::map(stats_list, function(x) max(x$dataset_id)))
  shifts <- c(0, max_ids[1:(length(max_ids)) - 1]) # Shift of IDs per dataset

  shift_dataset_id <- function(x, shift) {
    if(is.null(x)) {
      x
    } else {
      dplyr::mutate(x, dataset_id = dataset_id + shift)
    }
  }

  # Combines multiple data frame objects and then sorts by dataset_id
  bind_and_rearrange_df <- function(df_list) {
    dplyr::arrange(
      do.call(rbind, df_list),
      dataset_id
    )
  }

  # Apply the shifts of IDs to individual stats/diagnostics data frames
  stats_list <- purrr::map2(stats_list, shifts, shift_dataset_id)
  backend_diagnostics_list <- purrr::map2(backend_diagnostics_list, shifts, shift_dataset_id)
  default_diagnostics_list <- purrr::map2(default_diagnostics_list, shifts, shift_dataset_id)

  # Combine all the elements into a bigger object
  SBC_results(stats = bind_and_rearrange_df(stats_list),
              fits = do.call(c, fits_list),
              backend_diagnostics = bind_and_rearrange_df(backend_diagnostics_list),
              default_diagnostics = bind_and_rearrange_df(default_diagnostics_list),
              errors =  do.call(c, errors_list),
              messages = do.call(c, messages_list),
              warnings = do.call(c, warnings_list),
              outputs = do.call(c, outputs_list)
  )
}


#' @export
length.SBC_results <- function(x) {
  validate_SBC_results(x)
  length(x$fits)
}

#' Subset the results.
#'
#' @param indices integer indices or a binary vector of the same length as the number fits,
#' selecting which fits to keep.
#' @export
`[.SBC_results` <- function(x, indices) {
  validate_SBC_results(x)
  if(length(x) == 0 && length(indices) != 0) {
    stop("Cannot subset empty results with non-empty indices")
  }
  indices_to_keep <- (1:length(x))[indices]
  index_map <- 1:length(indices_to_keep)
  names(index_map) <- indices_to_keep

  subset_run_df <- function(df) {
    if(is.null(df)) {
      NULL
    }
    filtered <- dplyr::filter(df, dataset_id %in% indices_to_keep)
    remapped <- dplyr::mutate(filtered, dataset_id = index_map[as.character(dataset_id)])
    dplyr::arrange(remapped, dataset_id)
  }

  SBC_results(stats = subset_run_df(x$stats),
              fits = x$fits[indices],
              backend_diagnostics = subset_run_df(x$backend_diagnostics),
              default_diagnostics = subset_run_df(x$default_diagnostics),
              outputs = x$output[indices],
              messages = x$messages[indices],
              warnings = x$warnings[indices],
              errors = x$errors[indices])
}


#' Fit datasets and evaluate diagnostics and SBC metrics.
#'
#' Performs the main SBC routine given datasets and a backend.
#'
#' # Paralellization
#'
#' Parallel processing is supported via the `future` package, for most uses, it is most sensible
#'  to just call `plan(multisession)` once in your R session and  all
#'  cores your computer will be used. For more details refer to the documentation
#'  of the `future` package.
#'
#' # Thinning
#'
#' When using backends based on MCMC, there are two possible moments when
#' samples may need to be thinned. They can be thinned directly within the backend
#' and they may be thinned only to compute the ranks for SBC as specified by the
#' `thin_ranks` argument. The main reason those are separate is that computing the
#' ranks requires no or negligible autocorrelation while some autocorrelation
#' may be easily tolerated for summarising the fit results or assessing convergence.
#' In fact, thinning too aggressively in the backend may lead to overly noisy
#' estimates of posterior means, quantiles and the [posterior::rhat()] and
#' [posterior::ess_tail()] diagnostics. So for well-adapted Hamiltonian Monte-Carlo
#' chains (e.g. Stan-based backends), we recommend no thinning in the backend and
#' even value of `thin_ranks` between 6 and 10 is usually sufficient to remove
#' the residual autocorrelation. For a backend based on Metropolis-Hastings,
#' it might be sensible to thin quite aggressively already in the backend and
#' then have some additional thinning via `thin_ranks`.
#'
#' Backends that don't require thining should implement [SBC_backend_iid_samples()]
#' or [SBC_backend_default_thin_ranks()] to avoid thinning by default.
#'
#' @param datasets an object of class `SBC_datasets`
#' @param backend the model + sampling algorithm. The built-in backends can be constructed
#'   using [SBC_backend_cmdstan_sample()], [SBC_backend_cmdstan_variational()], [SBC_backend_rstan_sample()] and [SBC_backend_brms()].
#'   (more to come: issue 31, 38, 39). The backend is an S3 class supporting at least the [SBC_fit()],
#'   [SBC_fit_to_draws_matrix()] methods.
#' @param cores_per_fit how many cores should the backend be allowed to use for a single fit?
#'    Defaults to the maximum number that does not produce more parallel chains
#'    than you have cores. See [default_cores_per_fit()].
#' @param keep_fits boolean, when `FALSE` full fits are discarded from memory -
#'    reduces memory consumption and increases speed (when processing in parallel), but
#'    prevents you from inspecting the fits and using [recompute_statistics()].
#'    We recommend to set to `TRUE` in early phases of workflow, when you run just a few fits.
#'    Once the model is stable and you want to run a lot of iterations, we recommend setting
#'    to `FALSE` (even for quite a simple model, 1000 fits can easily exhaust 32GB of RAM).
#' @param thin_ranks how much thinning should be applied to posterior samples before computing
#'    ranks for SBC. Should be large enough to avoid any noticeable autocorrelation of the
#'    thinned samples.
#' @param chunk_size How many fits of `datasets` shall be processed in one batch
#'    by the same worker. Relevant only when using parallel processing.
#'    The larger the value, the smaller overhead there will be for parallel processing, but
#'    the work may be distributed less equally across workers. We recommend setting this high
#'    enough that a single batch takes at least several seconds, i.e. for small models,
#'    you can often reduce computation time noticeably by increasing this value.
#'    You can use `options(SBC.min_chunk_size = value)` to set a minimum chunk size globally.
#'    See documentation of `future.chunk.size` argument for [future.apply::future_lapply()] for more details.
#' @param cache_mode Type of caching of results, currently the only supported modes are
#'    `"none"` (do not cache) and `"results"` where the whole results object is stored
#'    and recomputed only when the hash of the backend or dataset changes.
#' @param cache_location The filesystem location of cache. For `cache_mode = "results"`
#'    this should be a name of a single file. If the file name does not end with
#'    `.rds`, this extension is appended.
#' @param globals A list of names of objects that are defined
#' in the global environment and need to present for the backend to work (
#' if they are not already available in package).
#' It is added to the `globals` argument to [future::future()], to make those
#' objects available on all workers.
#' @return An object of class [SBC_results()].
#' @export
compute_results <- function(datasets, backend,
                            cores_per_fit = default_cores_per_fit(length(datasets)),
                            keep_fits = TRUE,
                            thin_ranks = SBC_backend_default_thin_ranks(backend),
                            chunk_size = default_chunk_size(length(datasets)),
                            gen_quants = NULL,
                            cache_mode = "none",
                            cache_location = NULL,
                            globals = list()) {
  stopifnot(length(datasets) > 0)

  datasets <- validate_SBC_datasets(datasets)
  if(!is.null(gen_quants)) {
    gen_quants <- validate_generated_quantities(gen_quants)
  }

  ## Handle caching
  if(cache_mode == "results") {
    if(is.null(cache_location) || !dir.exists(dirname(cache_location))) {
      stop(SBC_error("SBC_invalid_argument_error",
                     "When using cache_mode == 'results', the cache_location argument must provide a filename in an existing directory"))
    }
    cache_basename <- basename(cache_location)
    if(!endsWith(cache_basename, ".rds")) {
      cache_location <- file.path(dirname(cache_location), paste0(cache_basename, ".rds"))
    }

    backend_hash <- SBC_backend_hash_for_cache(backend)
    data_hash <- rlang::hash(datasets)

    if(file.exists(cache_location)) {
      results_from_cache <- readRDS(cache_location)
      if(!is.list(results_from_cache) ||
         !all(
           c("result", "backend_hash", "data_hash", "thin_ranks", "gen_quants","keep_fits")
           %in% names(results_from_cache))) {
        warning("Cache file exists but is in invalid format. Will recompute.")
      } else if(results_from_cache$backend_hash != backend_hash) {
        message("Cache file exists but the backend hash differs. Will recompute.")
      } else if(results_from_cache$data_hash != data_hash) {
        message("Cache file exists but the datasets hash differs. Will recompute.")
      } else {
        result <- tryCatch(validate_SBC_results(results_from_cache$result),
                           error = function(e) { NULL })
        if(is.null(result)) {
          warning("Cache file contains invalid SBC_results object. Will recompute.")
        } else if(results_from_cache$thin_ranks != thin_ranks ||
                  !identical(results_from_cache$gen_quants, gen_quants))  {
          if(!results_from_cache$keep_fits) {
            message("Cache file exists, but was computed with different thin_ranks/gen_quants and keep_fits == FALSE. Will recompute.")
          } else {
            message(paste0("Results loaded from cache file '", cache_basename,
                           "' but it was computed with different thin_ranks/gen_quants.\n",
                           "Calling recompute_statistics."))
            return(recompute_statistics(old_results = result, datasets = datasets,
                                        thin_ranks = thin_ranks, gen_quants = gen_quants,
                                        backend = backend))
          }
        } else {
          message(paste0("Results loaded from cache file '", cache_basename, "'"))
          check_all_SBC_diagnostics(result)

          return(result)
        }
      }
    }
  } else if(cache_mode == "none") {
    if(!is.null(cache_location)) {
      warning("cache_location is provided, but cache_mode is set to 'none' - no caching will take place.")
    }
  } else {
    stop(SBC_error("SBC_invalid_argument_error", "Unrecognized cache mode"))
  }
  ## End of caching


  # Create combined data for computation
  params_and_generated_list <- list()
  for(i in 1:length(datasets)) {
    params_and_generated_list[[i]] <- list(
      parameters = posterior::subset_draws(datasets$parameters,
                                           draw = i),
      generated = datasets$generated[[i]]
    )
  }
  if(is.null(gen_quants)) {
    future.globals <- globals
  } else {
    future.globals <- c(globals, attr(gen_quants, "globals"))
  }

  results_raw <- future.apply::future_lapply(
    params_and_generated_list, SBC:::compute_results_single,
    backend = backend, cores = cores_per_fit,
    keep_fit = keep_fits, thin_ranks = thin_ranks,
    gen_quants = gen_quants,
    future.seed = TRUE,
    future.globals = future.globals,
    future.chunk.size = chunk_size)

  # Combine, check and summarise
  fits <- rep(list(NULL), length(datasets))
  outputs <- rep(list(NULL), length(datasets))
  messages <- rep(list(NULL), length(datasets))
  warnings <- rep(list(NULL), length(datasets))
  errors <- rep(list(NULL), length(datasets))
  stats_list <- list()
  backend_diagnostics_list <- list()
  n_errors <- 0
  max_errors_to_show <- 5
  for(i in 1:length(datasets)) {
    if(!is.null(results_raw[[i]]$fit)) {
      fits[[i]] <- results_raw[[i]]$fit
    }
    if(is.null(results_raw[[i]]$error)) {
      stats_list[[i]] <- results_raw[[i]]$stats
      stats_list[[i]]$dataset_id <- i
      stats_list[[i]] <- dplyr::select(stats_list[[i]], dataset_id, tidyselect::everything())
      backend_diagnostics_list[[i]] <- results_raw[[i]]$backend_diagnostics
      if(!is.null(results_raw[[i]]$backend_diagnostics)){
        backend_diagnostics_list[[i]]$dataset_id <- i
        backend_diagnostics_list[[i]] <- dplyr::select(backend_diagnostics_list[[i]], dataset_id, tidyselect::everything())
      }
    }
    else {
      if(n_errors < max_errors_to_show) {
        if(is.null(results_raw[[i]]$fit)) {
          message("Dataset ", i, " resulted in error when fitting.\n")
          message(results_raw[[i]]$error, "\n")
          if(!is.null(results_raw[[i]]$warnings)) {
            message(" --- Warnings for fit ", i, " ----")
            message(paste0(results_raw[[i]]$warnings, collapse = "\n"))
          }
          if(!is.null(results_raw[[i]]$messages)) {
            message(" --- Messages for fit ", i, " ----")
            message(paste0(results_raw[[i]]$messages, collapse = "\n"))
          }
          if(is.null(results_raw[[i]]$output)) {
            message(" --- Nothing in stdout ---")
          } else {
            message(" ---- Model output ----")
            cat(paste0(results_raw[[i]]$output, collapse = "\n"))
          }
          message("\n ---- End of output for dataset ", i, " -----")
        } else {
          message("Dataset ", i, " resulted in error when post-processing the fit.\n",
                  "Calling `recompute_statistics` after you've found and fixed the problem could ",
                  "let you move further without refitting")
          message(results_raw[[i]]$error, "\n")
        }

      } else if(n_errors == max_errors_to_show) {
        message("Too many datasets produced errors. Further error messages not shown.\n")
      }
      n_errors <- n_errors + 1
      errors[[i]] <- results_raw[[i]]$error
    }
    if(!is.null(results_raw[[i]]$output)) {
      outputs[[i]] <- results_raw[[i]]$output
    }
    if(!is.null(results_raw[[i]]$messages)) {
      messages[[i]] <- results_raw[[i]]$messages
    }
    if(!is.null(results_raw[[i]]$warnings)) {
      warnings[[i]] <- results_raw[[i]]$warnings
    }
  }

  if(n_errors == length(datasets)) {
    warning("All datasets produced error when fitting")
  } else if(n_errors > 0) {
    warning("Total of ", n_errors, " datasets produced errors.")
  }

  stats <- do.call(rbind, stats_list)
  backend_diagnostics <- do.call(rbind, backend_diagnostics_list)

  if(!is.null(stats)) {
    check_stats(stats, datasets, thin_ranks)
  } else {
    # Return dummy stats that let the rest of the code work.
    stats <- data.frame(dataset_id = integer(0), rhat = numeric(0), ess_bulk = numeric(0),
                        ess_tail = numeric(0),
                        rank = integer(0), simulated_value = numeric(0), max_rank = integer(0))
  }

  default_diagnostics <-  tryCatch(
    { compute_default_diagnostics(stats) },
    error = function(e) { warning("Error when computing param diagnostics. ", e); NULL })


  res <- SBC_results(stats = stats, fits = fits, outputs = outputs,
                     messages = messages,
                     warnings = warnings,
                     backend_diagnostics = backend_diagnostics,
                     default_diagnostics = default_diagnostics,
                     errors = errors)

  if(cache_mode == "results") {
    results_for_cache <- list(result = res, backend_hash = backend_hash,
                              data_hash = data_hash, thin_ranks = thin_ranks,
                              gen_quants = gen_quants, keep_fits = keep_fits)
    tryCatch(saveRDS(results_for_cache, file = cache_location),
             error = function(e) { warning("Error when saving cache file: ", e) })
  }

  check_all_SBC_diagnostics(res)

  res
}

#' Determines the default cores per single fit.
#'
#' When parallel processing is disabled, this just returns the number of available
#' cores. Otherwise, it chooses the largest integer that keeps
#' `cores_per_fit * (n_fits / chunk_size) <= total_cores`, i.e. it avoids
#' running so many chains in parallel that there will be more chains than cores.
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

#' Determines the default chunk size.
#'
#' By default will make every worker process a single chunk.
#' You can set the `options(SBC.min_chunk_size = value)` to enforce a minimum
#' chunk size globally (chunk size can still be larger if you have substantially more
#' fits to run than workers.
#' @export
default_chunk_size <- function(n_fits, n_workers = future::nbrOfWorkers()) {
  guess <- if(is.infinite(n_workers)) {
    1
  } else {
    n_fits / n_workers
  }
  max(guess, getOption("SBC.min_chunk_size", 1))
}

# Capturing output.
# Based on https://www.r-bloggers.com/2020/10/capture-message-warnings-and-errors-from-a-r-function/
capture_all_outputs <- function(expr) {
  logs <- list(message = list(), warning = list())
  add_log <- function(type, message) {
    new_l <- logs
    new_l[[type]][[length(new_l[[type]]) + 1]]  <- message
    logs <<- new_l
  }
  output <- capture.output({
    res <- withCallingHandlers(
      expr,
      warning=function(w) {
        add_log("warning", conditionMessage(w))
        invokeRestart("muffleWarning")
      }, message = function(m) {
        add_log("message", conditionMessage(m))
        invokeRestart("muffleMessage")
      })
  }, type = "output")
  list(result = res, messages = do.call(c, logs$message), warnings = do.call(c, logs$warning), output = output)
}


compute_results_single <- function(params_and_generated, backend, cores,
                                   keep_fit, thin_ranks, gen_quants) {

  parameters <- params_and_generated$parameters
  generated <- params_and_generated$generated

  # Note: explicitly referencing functions from the SBC package is needed
  # here as the function might be run in a separate R session that does not
  # have SBC imported.
  result_with_output <- SBC:::capture_all_outputs({
    res <- tryCatch({
      fit <- SBC_fit(backend, generated, cores = cores)
      c(list(fit = fit, error = NULL))
    }, error = function(e) { list(fit = NULL, error = e) })
  })

  res <- result_with_output$result

  res$output <- result_with_output$output
  res$messages <- result_with_output$messages
  res$warnings <- result_with_output$warnings

  if(is.null(res$error)) {
    error_stats <-  SBC:::capture_all_outputs({
      tryCatch( {
        res$stats <- SBC::statistics_from_single_fit(
          res$fit, parameters = parameters, thin_ranks = thin_ranks,
          generated = generated, gen_quants = gen_quants,
          backend = backend)

        res$backend_diagnostics <- SBC::SBC_fit_to_diagnostics(
          fit, res$output, res$messages, res$warnings)
        NULL
      }, error = identity)
    })

    if(!is.null(error_stats$result)) {
      res$error <- error_stats$result
    }

    if(!is.null(error_stats$output) && length(error_stats$output > 0)) {
      res$output <- c(res$output, "\n== Output from computing statistics ==\n", error_stats$output)
    }
    if(!is.null(error_stats$messages) && length(error_stats$messages > 0)) {
      res$messages <- c(res$message, "== Messages from computing statistics ==", error_stats$messages)
    }
    if(!is.null(error_stats$warnings) && length(error_stats$warnings > 0)) {
      res$warnings <- c(res$warnings, "== Warnings from computing statistics ==", error_stats$warnings)
    }
  } else {
    res$stats <- NULL
    res$backend_diagnostics <- NULL
  }

  if(!keep_fit) {
    res$fit <- NULL
  }

  res
}

#' Recompute SBC statistics given a single fit.
#'
#' Potentially useful for doing some advanced stuff, but should not
#' be used in regular workflow. Use [recompute_statistics()] to update
#' an `[SBC_results]` objects with different `thin_ranks` or other settings.
#'
#' @export
#' @seealso [recompute_statistics()]
statistics_from_single_fit <- function(fit, parameters, generated,
                                       thin_ranks, gen_quants,
                                       backend) {

  fit_matrix <- SBC_fit_to_draws_matrix(fit)

  if(!is.null(gen_quants)){
    gen_quants <- validate_generated_quantities(gen_quants)
    gq_fit <- compute_gen_quants(fit_matrix, generated, gen_quants)
    fit_matrix <- posterior::bind_draws(fit_matrix, gq_fit, along = "variable")

    gq_parameter <- compute_gen_quants(parameters, generated, gen_quants)
    parameters <- posterior::bind_draws(parameters, gq_parameter, along = "variable")
  }

  shared_pars <- intersect(posterior::variables(parameters),
                           posterior::variables(fit_matrix))


  # Make sure the order of parameters matches
  parameters <- posterior::subset_draws(parameters, variable = shared_pars)


  fit_matrix <- posterior::subset_draws(fit_matrix, variable = shared_pars)

  fit_thinned <- posterior::thin_draws(fit_matrix, thin_ranks)


  stats <- posterior::summarise_draws(fit_matrix)

  if(SBC_backend_iid_samples(backend)) {
    ## iid samples have the bestest diagnostics by construction
    stats$rhat <- 1
    stats$ess_bulk <- posterior::ndraws(fit_matrix)
    stats$ess_tail <- posterior::ndraws(fit_matrix)
  }

  stats <- dplyr::rename(stats, parameter = variable)
  stats$simulated_value <- as.numeric(parameters)

  ranks <- calculate_ranks_draws_matrix(parameters, fit_thinned)
  if(!identical(stats$parameter, names(ranks))) {
    stop("A naming conflict")
  }
  stats$rank <- ranks
  stats$max_rank <- attr(ranks, "max_rank")
  stats$z_score <- (stats$simulated_value - stats$mean) / stats$sd

  stats <- dplyr::select(
    stats, parameter, simulated_value, rank, z_score, tidyselect::everything())

  stats
}

# check that the computed stats data frame hs problems
check_stats <- function(stats, datasets, thin_ranks) {
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

  all_pars <- dplyr::summarise(
    dplyr::group_by(stats, dataset_id),
    all_pars = paste0(parameter, collapse = ","), .groups = "drop")
  if(length(unique(all_pars$all_pars)) != 1) {
    warning("Not all fits share the same parameters")
  }

  missing_pars <- setdiff(posterior::variables(datasets$parameters), stats$parameter)
  if(length(missing_pars) > 0) {
    warning("Some parameters missing in fits: ", paste0(missing_pars, collapse = ", "))

  }
}

#' Create a definition of generated quantities evaluated in R.
#'
#' When the expression contains non-library functions/objects, and parallel processing
#' is enabled, those must be
#' named in the `.globals` parameter (hopefully we'll be able to detect those
#' automatically in the future). Note that [recompute_statistics()] currently
#' does not use parallel processing, so `.globals` don't need to be set.
#'
#' @param ... named expressions representing the quantitites
#' @param .globals A list of names of objects that are defined
#' in the global environment and need to present for the gen. quants. to evaluate.
#' It is added to the `globals` argument to [future::future()], to make those
#' objects available on all workers.
#' @export
generated_quantities <- function(..., .globals = list()) {
  structure(rlang::enquos(..., .named = TRUE),
            class = "SBC_generated_quantities",
            globals = .globals
            )
}

#' @export
validate_generated_quantities <- function(x) {
  stopifnot(inherits(x, "SBC_generated_quantities"))
  invisible(x)
}

#'@export
compute_gen_quants <- function(draws, generated, gen_quants) {
  gen_quants <- validate_generated_quantities(gen_quants)
  draws_rv <- posterior::as_draws_rvars(draws)

  draws_env <- list2env(draws_rv)
  generated_env <- list2env(generated, parent = draws_env)

  data_mask <- rlang::new_data_mask(bottom = generated_env, top = draws_env)

  eval_func <- function(gq) {
    # Wrap the expression in `rdo` which will mostly do what we need
    # all the tricks are just to have the correct environment when we need it
    wrapped_gq <- rlang::new_quosure(rlang::expr(posterior::rdo(!!rlang::get_expr(gq))), rlang::get_env(gq))
    rlang::eval_tidy(wrapped_gq, data = data_mask)
  }
  rvars <- lapply(gen_quants, FUN = eval_func)
  do.call(posterior::draws_rvars, rvars)
}

#' Recompute SBC statistics without refitting models.
#'
#' Useful for example to recompute SBC ranks with a different choice of `thin_ranks`
#' or added generated quantities.
#' @return An S3 object of class `SBC_results` with updated `$stats` and `$default_diagnostics` fields.
#' @param backend backend used to fit the results. Used to pull various defaults
#'   and other setting influencing the computation of statistics.
#' @export
recompute_statistics <- function(old_results, datasets, backend,
                                 thin_ranks = SBC_backend_default_thin_ranks(backend),
                                 gen_quants = NULL) {
  validate_SBC_results(old_results)
  validate_SBC_datasets(datasets)

  if(length(old_results) != length(datasets)) {
    stop("The number of fits in old_results does not match the number of datasets")
  }

  new_results <- old_results
  missing_fits <- purrr::map_lgl(old_results$fits, is.null)
  if(all(missing_fits)) {
    stop("No raw fits preserved, cannot recompute. ",
         "Either all datasets produced errors or the results were computed with keep_fits = FALSE")
  } else if(any(missing_fits)) {
    warning("Some raw fits not available. Those fits will be ignored when recomputing statistics")
  }

  new_stats_list <- list()
  for(i in 1:length(old_results)) {
    if(!is.null(old_results$fits[[i]])) {
      parameters <- posterior::subset_draws(datasets$parameters, draw = i)
      new_stats_list[[i]] <- statistics_from_single_fit(old_results$fits[[i]],
                                                        parameters = parameters,
                                                        generated = datasets$generated[[i]],
                                                        thin_ranks = thin_ranks,
                                                        gen_quants = gen_quants,
                                                        backend = backend)
      new_stats_list[[i]]$dataset_id <- i
      new_stats_list[[i]] <- dplyr::select(new_stats_list[[i]], dataset_id, tidyselect::everything())

    }
  }

  new_stats <- do.call(rbind, new_stats_list)
  check_stats(new_stats, datasets, thin_ranks)

  new_results$stats <- new_stats

  new_results$default_diagnostics <-  tryCatch(
    { compute_default_diagnostics(new_stats) },
    error = function(e) { warning("Error when computing param diagnostics. ", e); NULL })


  check_all_SBC_diagnostics(new_results)

  new_results

}

#' Discrete uniform distribution allowing for varying lower and upper bounds.
#'
#' Internal, should not be exported.
#' Based on https://stats.stackexchange.com/a/3939/73129
#' @keywords internal
rdunif <- function(n, a, b) {
  ceiling(runif(n, min = a - 1, max= b))
}

#' Calculate ranks given parameter values within a posterior distribution.
#'
#' When there are ties (e.g. for discrete parameters), the rank is currently drawn stochastically
#' among the ties.
#' @param params a vector of values to check
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
SBC_diagnostic_messages <- function(message_df) {
  if(!inherits(message_df, "SBC_diagnostic_messages")) {
    class(message_df) <- c("SBC_diagnostic_messages", class(message_df))
  }
  validate_diagnostic_messages(message_df)
}

validate_diagnostic_messages <- function(x) {
  stopifnot(is.data.frame(x))
  stopifnot(inherits(x, "SBC_diagnostic_messages"))
  if(!identical(sort(names(x)), c("message", "ok"))) {
    stop("Diagnostic messages have to have columns 'message' and 'ok' and no other")
  }

  x
}

print.SBC_diagnostic_messages <- function(x, include_ok = TRUE, print_func = cat) {
  x <- validate_diagnostic_messages(x)
  if(!include_ok) {
    x <- dplyr::filter(x, !ok)
  }

  if(nrow(x) > 0) {
    for(i in 1:nrow(x)) {
      print_func(paste0(" - ", x$message[i], "\n"))
    }
  }
}

#' Get diagnostic messages for `SBC_results` or other objects.
#'
#' @export
#' @return An object of class `SBC_diagnostic_messages`, inheriting a data.frame.
get_diagnostic_messages <- function(x) {
  UseMethod("get_diagnostic_messages")
}


#' Check diagnostics and issue warnings when those fail.
#'
#' @rdname check_all_SBC_diagnostics
#' @export
#' @return TRUE if all checks are OK, FALSE otherwise.
check_all_SBC_diagnostics <- function(x) {
  UseMethod("check_all_SBC_diagnostics")
}

#' @export
#' @rdname check_all_SBC_diagnostics
check_all_SBC_diagnostics.default <- function(x) {
  if(!is.null(x)) {
    msg <- get_diagnostic_messages(x)
    print(msg, include_ok = FALSE, print_func = function(m) { message(m, appendLF = FALSE) })
    invisible(all(msg$ok))
  } else {
    invisible(TRUE)
  }
}

#' @rdname check_all_SBC_diagnostics
#' @export
check_all_SBC_diagnostics.SBC_results <- function(x) {
  res <- NextMethod()
  if(!res) {
    message("Not all diagnostics are OK.\nYou can learn more by inspecting $default_diagnostics, ",
            "$backend_diagnostics \nand/or investigating $outputs/$messages/$warnings for detailed output from the backend.")
  }
  res
}

#' @export
get_diagnostic_messages.SBC_results <- function(x) {
  get_diagnostic_messages(summary(x))
}



#' @export
summary.SBC_results <- function(x) {
  summ <- list(
    n_fits = length(x$fits),
    n_errors = sum(!purrr::map_lgl(x$errors, is.null)),
    n_warnings = sum(purrr::map_lgl(x$messages, ~ !is.null(.x) && any(x$type == "warning"))),
    n_high_rhat = sum(is.na(x$default_diagnostics$max_rhat) | x$default_diagnostics$max_rhat > 1.01),
    max_max_rhat = max(c(-Inf, x$default_diagnostics$max_rhat)),
    n_low_ess_to_rank = sum(is.na(x$default_diagnostics$min_ess_to_rank) | x$default_diagnostics$min_ess_to_rank < 0.5),
    min_min_ess_bulk = min(c(Inf, x$default_diagnostics$min_ess_bulk)),
    min_min_ess_tail = min(c(Inf, x$default_diagnostics$min_ess_tail))
  )
  if(!is.null(x$backend_diagnostics)) {
    summ$backend_diagnostics <- summary(x$backend_diagnostics)
  } else {
    summ$backend_diagnostics <- NULL
  }
  structure(
    summ,
    class = "SBC_results_summary"
  )
}

#' @export
get_diagnostic_messages.SBC_results_summary <- function(x) {

  message_list <- list()
  i <- 1
  if(x$n_errors > 0) {
    msg <- paste0(x$n_errors, " (", round(100 * x$n_errors / x$n_fits), "%) fits resulted in an error.")
    message_list[[i]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had errors.")
  }
  i <- i + 1

  if(x$n_warnings > 0) {
    msg <- paste0(x$n_warnings, " (", round(100 * x$n_warnings / x$n_fits), "%) fits gave warnings. Inspect $messages to see them.")
    message_list[[i]] <- data.frame(ok = TRUE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "No fits gave warnings.")
  }
  i <- i + 1

  if(x$n_high_rhat > 0) {
    msg <- paste0(x$n_high_rhat, " (", round(100 * x$n_high_rhat / x$n_fits), "%) fits had at least one Rhat > 1.01. ",
                  "Largest Rhat was ",round(x$max_max_rhat, 3), ".")
    message_list[[i]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had Rhat > 1.01.")
  }
  i <- i + 1

  if(x$n_low_ess_to_rank > 0) {
    msg <- paste0(x$n_low_ess_to_rank, " (", round(100 * x$n_low_ess_to_rank / x$n_fits), "%) fits had tail ESS undefined or less than ",
                  "half of the maximum rank, potentially skewing \nthe rank statistics. The lowest tail ESS was ", round(x$min_min_ess_tail),
                  ".\n If the fits look good otherwise, increasing `thin_ranks` (via recompute_statistics) \nor number of posterior samples (by refitting) might help.")
    message_list[[i]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "All fits had tail ESS > half of the maximum rank.")
  }
  i <- i + 1

  message_list[[i]] <- data.frame(ok = TRUE, message = paste0("The lowest bulk ESS was ", round(x$min_min_ess_bulk)))
  i <-  i + 1

  if(!is.null(x$backend_diagnostics)) {
    message_list[[i]] <- get_diagnostic_messages(x$backend_diagnostics)
    i <- i + 1
  }

  SBC_diagnostic_messages(do.call(rbind, message_list))
}

#' @export
print.SBC_results_summary <- function(x) {
  cat("SBC_results with", x$n_fits, "total fits.\n")

  msg <- get_diagnostic_messages(x)
  print(msg)

  if(!all(msg$ok)) {
    message("Not all diagnostics are OK.\nYou can learn more by inspecting $default_diagnostics, ",
            "$backend_diagnostics \nand/or investigating $outputs/$messages/$warnings for detailed output from the backend.")
  }


  invisible(x)
}
