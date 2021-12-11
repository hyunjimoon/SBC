#' Use backend to fit a model to data.
#'
#' S3 generic, needs to be implemented by all backends.
#' All implementations have to return an object for which you can safely
#' call [SBC_fit_to_draws_matrix()] and get some draws.
#' If that's not possible an error should be raised.
#' @export
SBC_fit <- function(backend, generated, cores) {
  UseMethod("SBC_fit")
}

#' @export
SBC_fit_to_draws_matrix <- function(fit) {
  UseMethod("SBC_fit_to_draws_matrix")
}

#' @export
SBC_fit_to_draws_matrix.default <- function(fit) {
  posterior::as_draws_matrix(fit)
}

#' S3 generic to get backend-specific diagnostics.
#'
#' The diagnostics object has to be a `data.frame` but may
#' inherit additional classes - in particular it may be useful
#' for the returning object to implement [get_diagnostic_messages()].
#'
#' @param fit The fit returned by `SBC_fit`
#' @param fit_output a character string capturing what the backend wrote to stdout
#' @param fit_messages a character vector of messages the backend raised
#' @param fit_warnings a character vector of warnings the backend raised
#' @return an single row `data.frame` that includes diagnostics or NULL, if no diagnostics available.
#' @export
SBC_fit_to_diagnostics <- function(fit, fit_output, fit_messages, fit_warnings) {
  UseMethod("SBC_fit_to_diagnostics")
}

#' @export
SBC_fit_to_diagnostics.default <- function(fit, fit_output, fit_messages, fit_warnings) {
  NULL
}

#' Get hash used to identify cached results.
#'
#' S3 generic that allows backends to override how a hash is computed. By default `rlang::hash()`
#' is used.
#'
#' @export
SBC_backend_hash_for_cache <- function(backend) {
  UseMethod("SBC_backend_hash_for_cache")
}

#' @export
SBC_backend_hash_for_cache.default <- function(backend) {
  rlang::hash(backend)
}

#' S3 generic to let backends signal that they produced independent samples.
#'
#' Most backends (e.g. those based on variatns of MCMC) don't produce
#' independent samples and thus diagnostics like Rhat and ESS are important
#' and samples may need thinning. Backends that already produce independent
#' samples (e.g. ADVI/optimizing) can implement this method to return `TRUE`
#' to signal this is the case. The default implementation returns `FALSE`.
#' @param backend to check
#' @export
SBC_backend_iid_samples <- function(backend) {
  UseMethod("SBC_backend_iid_samples")
}

#' @rdname SBC_backend_iid_samples
#' @export
SBC_backend_iid_samples.default <- function(backend) {
  FALSE
}

#' S3 generic to get backend-specific default thinning for rank computation.
#'
#' The default implementation plays it relatively safe and returns 10, unless
#' [SBC_backend_iid_samples()] returns `TRUE` in which case it returns 1.
#'
#' @export
SBC_backend_default_thin_ranks <- function(backend) {
  UseMethod("SBC_backend_default_thin_ranks")
}

#' @rdname SBC_backend_default_thin_ranks
#' @export
SBC_backend_default_thin_ranks.default <- function(backend) {
  if(SBC_backend_iid_samples(backend)) {
    1
  } else {
    10
  }
}



#' SBC backend using the `sampling` method from `rstan`.
#'
#' @param model a `stanmodel` object (created via `rstan::stan_model`)
#' @param ... other arguments passed to `sampling` (number of iterations, ...).
#'   Arguments `data` and `cores` cannot be set this way as they need to be
#'   controlled by the package.
#' @export
SBC_backend_rstan_sample <- function(model, ...) {
  stopifnot(inherits(model, "stanmodel"))
  args <- list(...)
  unacceptable_params <- c("data", "cores")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }

  if(is.null(args$open_progress)) {
    args$open_progress <- FALSE
  }
  structure(list(model = model, args = args), class = "SBC_backend_rstan_sample")
}


#' @export
SBC_fit.SBC_backend_rstan_sample <- function(backend, generated, cores) {
  do.call(rstan::sampling,
          combine_args(list(object = backend$model,
                 data = generated,
                 ## TODO: Forcing a single core until we can capture output with multiple cores
                 ## https://discourse.mc-stan.org/t/capturing-warnings-rejects-from-rstan-with-multiple-cores/23976
                 cores = 1),
            backend$args
            ))
}

#' @export
SBC_fit_to_diagnostics.stanfit <- function(fit, fit_output, fit_messages, fit_warnings) {
  res <- data.frame(
    max_chain_time = max(rowSums(rstan::get_elapsed_time(fit))),
    n_divergent = rstan::get_num_divergent(fit),
    n_max_treedepth = rstan::get_num_max_treedepth(fit),
    min_bfmi = min(rstan::get_bfmi(fit)),
    n_rejects = sum(grepl("reject", fit_messages)) + sum(grepl("reject", fit_warnings))
  )

  class(res) <- c("SBC_nuts_diagnostics", class(res))
  res
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_rstan_sample <- function(backend) {
  rlang::hash(list(model = backend$model@model_code, args = backend$args))
}

#' @export
summary.SBC_nuts_diagnostics <- function(diagnostics) {
  summ <- list(
    n_fits = nrow(diagnostics),
    max_chain_time = max(diagnostics$max_chain_time),
    has_divergent = sum(diagnostics$n_divergent > 0),
    max_divergent = max(diagnostics$n_divergent),
    has_treedepth = sum(diagnostics$n_max_treedepth > 0),
    max_max_treedepth = max(diagnostics$n_max_treedepth),
    has_rejects = sum(diagnostics$n_rejects > 0),
    max_rejects = max(diagnostics$n_rejects)
  )

  if(!is.null(diagnostics$min_bfmi)) {
    summ$has_low_bfmi = sum(diagnostics$min_bfmi < 0.2)
  }

  if(!is.null(diagnostics$n_failed_chains)) {
    if(any(is.na(diagnostics$n_failed_chains))) {
      problematic_fit_ids <- paste0(which(is.na(diagnostics$n_failed_chains)), collapse = ", ")
      warning("Fits for datasets ", problematic_fit_ids, " had NA for n_failed_chains.")
    }
    summ$has_failed_chains = sum(is.na(diagnostics$n_failed_chains) | diagnostics$n_failed_chains > 0)
  }

  structure(summ, class = "SBC_nuts_diagnostics_summary")
}


get_diagnostic_messages.SBC_nuts_diagnostics <- function(x) {
  get_diagnostic_messages(summary(x))
}

# Not currently used. Discussion at https://discourse.mc-stan.org/t/summarising-rhat-values-over-multiple-variables-fits/23957/
# was helpful and need to figure out a better way to handle this.
#
# Use the Generalized extreme value distribution
# to get a quantile of maximum of `n_vars` random values
# distributed as N(1, 0.005).
# Following https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution#Example_for_Normally_distributed_variables
# The approximation looks good for n_vars >= 10
# for smaller, we just plug in a log function with appropriate scale
# See https://math.stackexchange.com/questions/89030/expectation-of-the-maximum-of-gaussian-random-variables
# for a discussion on why log
get_expected_max_rhat <- function(n_vars, prob = 0.99, approx_sd = 0.005) {
  stopifnot(is.numeric(n_vars))
  stopifnot(all(n_vars >= 1))

  # Maximum of n_vars standardized normals
  gumbel_approx <- function(n) {
    # Gumbel params
    mu_n <- qnorm(1 - 1/n)
    sigma_n <- qnorm(1 - (1 / n) * exp(-1)) - mu_n

    # Inverse CDF of gumbel with xi = 0
    mu_n - (sigma_n * log(-log(prob)))
  }


  linear_bound <- 10
  approx_at_bound <- gumbel_approx(linear_bound)
  value_at_1 <- qnorm(prob)
  linear_scale <- (approx_at_bound - value_at_1 ) / log(linear_bound)

  std_val_max <- dplyr::if_else(n_vars < linear_bound,
    value_at_1 + linear_scale * log(n_vars),
    gumbel_approx(n_vars)
  )
  1 + std_val_max * approx_sd
}

get_diagnostic_messages.SBC_nuts_diagnostics_summary <- function(x) {
  message_list <- list()
  i <- 1
  if(!is.null(x$has_failed_chains)) {
    if(x$has_failed_chains > 0) {
      msg <- paste0(x$has_failed_chains, " (", round(100 * x$has_failed_chains / x$n_fits), "%) fits had some failed chains.")
      message_list[[i]] <- data.frame(ok = FALSE, message = msg)
    } else {
      message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had failed chains.")
    }
    i <- i + 1
  }

  if(x$has_divergent > 0) {
    msg <- paste0(x$has_divergent, " (", round(100 * x$has_divergent / x$n_fits),
                  "%) fits had divergent transitions. Maximum number of divergences was ",
                  x$max_divergent, ".")
    message_list[[i]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had divergent transitions.")
  }
  i <- i + 1

  if(x$has_treedepth > 0) {
    msg <- paste0(x$has_treedepth, " (", round(100 * x$has_treedepth / x$n_fits),
                  "%) fits had iterations that saturated max treedepth. Maximum number of max treedepth was ",
                  x$max_max_treedepth, ".")
    message_list[[i]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had iterations that saturated max treedepth.")
  }
  i <- i + 1

  if(!is.null(x$has_low_bfmi)) {
    if(x$has_low_bfmi > 0) {
      msg <- paste0(x$has_low_bfmi, " (", round(100 * x$has_low_bfmi / x$n_fits), "%) fits had low BFMI.")
      message_list[[i]] <- data.frame(ok = FALSE, message = msg)
    } else {
      message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had low BFMI.")
    }
    i <- i + 1
  }

  if(x$has_rejects > 0) {
    msg <- paste0(x$has_rejects, " (", round(100 * x$has_rejects / x$n_fits), "%) fits had some steps ",
                  "rejected. Maximum number of rejections was ", x$max_rejects, ".")
    message_list[[i]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had steps rejected.")
  }
  i <- i + 1

  message_list[[i]] <- data.frame(ok = TRUE, message = paste0("Maximum time per chain was ", x$max_chain_time, " sec."))
  i <- i + 1

  SBC_diagnostic_messages(do.call(rbind, message_list))
}

#' @export
print.SBC_nuts_diagnostics_summary <- function(x) {
  msg <- get_diagnostic_messages(x)
  print(msg)
  invisible(x)
}

#' Backend based on sampling via `cmdstanr`.
#'
#' @param model an object of class `CmdStanModel` (as created by `cmdstanr::cmdstan_model`)
#' @param ... other arguments passed to the `$sample()` method of the model. The `data` and
#'   `parallel_chains` arguments cannot be set this way as they need to be controlled by the SBC
#'   package.
#' @export
SBC_backend_cmdstan_sample <- function(model, ...) {
  require_cmdstanr_version("cmdstan backend")

  stopifnot(inherits(model, "CmdStanModel"))
  if(length(model$exe_file()) == 0) {
    stop("The model has to be already compiled, call $compile() first.")
  }
  args <- list(...)
  unacceptable_params <- c("data", "parallel_chains ", "cores", "num_cores")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }
  structure(list(model = model, args = args), class = "SBC_backend_cmdstan_sample")
}

#' @export
SBC_fit.SBC_backend_cmdstan_sample <- function(backend, generated, cores) {
  fit <- do.call(backend$model$sample,
          combine_args(backend$args,
            list(
              data = generated,
              parallel_chains = cores)))

  if(all(fit$return_codes() != 0)) {
    stop("No chains finished succesfully")
  }

  fit
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_cmdstan_sample <- function(backend) {
  rlang::hash(list(model = backend$model$code(), args = backend$args))
}


#' @export
SBC_fit_to_draws_matrix.CmdStanMCMC <- function(fit) {
  fit$draws(format = "draws_matrix")
}


#' @export
SBC_fit_to_diagnostics.CmdStanMCMC <- function(fit, fit_output, fit_messages, fit_warnings) {
  res <- data.frame(
    max_chain_time = max(fit$time()$chains[,"total"]),
    n_failed_chains = fit$num_chains() - sum(fit$return_codes() == 0),
    n_divergent = sum(fit$sampler_diagnostics()[, , "divergent__"]),
    n_max_treedepth =  sum(fit$sampler_diagnostics()[, , "treedepth__"] == fit$metadata()$max_treedepth),
    n_rejects = sum(grepl("reject", fit_messages)) + sum(grepl("reject", fit_warnings))
    #
  ) # TODO: add min_bfmi once https://github.com/stan-dev/cmdstanr/pull/500/ is merged
  class(res) <- c("SBC_nuts_diagnostics", class(res))
  res
}

#' Backend based on variational approximation via `cmdstanr`.
#'
#' @param model an object of class `CmdStanModel` (as created by `cmdstanr::cmdstan_model`)
#' @param ... other arguments passed to the `$variational()` method of the model.
#' The `data` argument cannot be set this way as they need to be controlled by the SBC
#'   package.
#' @export
SBC_backend_cmdstan_variational <- function(model, ...) {
  require_cmdstanr_version("cmdstan backend")

  stopifnot(inherits(model, "CmdStanModel"))
  if(length(model$exe_file()) == 0) {
    stop("The model has to be already compiled, call $compile() first.")
  }
  args <- list(...)
  unacceptable_params <- c("data")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }

  structure(list(model = model, args = args), class = "SBC_backend_cmdstan_variational")
}

#' @export
SBC_fit.SBC_backend_cmdstan_variational <- function(backend, generated, cores) {
  fit <- do.call(backend$model$variational,
                 combine_args(backend$args,
                              list(
                                data = generated)))

  if(all(fit$return_codes() != 0)) {
    stop("Variational inference did not finish succesfully")
  }

  fit
}


#' @export
SBC_backend_hash_for_cache.SBC_backend_cmdstan_variational <- function(backend) {
  rlang::hash(list(model = backend$model$code(), args = backend$args))
}


#' @export
SBC_fit_to_draws_matrix.CmdStanVB <- function(fit) {
  fit$draws(format = "draws_matrix")
}

#' @export
SBC_backend_iid_samples.SBC_backend_cmdstan_variational <- function(backend) {
  TRUE
}


#' @export
SBC_fit_to_diagnostics.CmdStanVB <- function(fit, fit_output, fit_messages, fit_warnings) {
  res <- data.frame(
    elbo_converged = any(grepl("ELBO CONVERGED", fit_output)),
    n_rejects = sum(grepl("reject", fit_messages)) + sum(grepl("reject", fit_warnings)),
    time = fit$time()$total
  )

  class(res) <- c("SBC_ADVI_diagnostics", class(res))
  res
}

#' @export
summary.SBC_ADVI_diagnostics <- function(x) {
  summ <- list(
    n_fits = nrow(x),
    n_elbo_not_converged = sum(!x$elbo_converged),
    has_rejects = sum(x$n_rejects > 0),
    max_rejects = max(x$n_rejects),
    max_time = max(x$time)

  )

  structure(summ, class = "SBC_ADVI_diagnostics_summary")
}

#' @export
get_diagnostic_messages.SBC_ADVI_diagnostics <- function(x) {
  get_diagnostic_messages(summary(x))
}


#' @export
get_diagnostic_messages.SBC_ADVI_diagnostics_summary <- function(x) {
  message_list <- list()

  if(x$n_elbo_not_converged == 0) {
    message_list[[1]] <-
      data.frame(ok = TRUE, message = "All fits converged.")
  } else {
    message_list[[1]] <-
      data.frame(ok = FALSE,
                 message = paste0(
                   x$n_elbo_not_converged, " (", round(100 * x$n_elbo_not_converged / x$n_fits),
                   "%) of fits did not converge."))
  }


  if(x$has_rejects > 0) {
    msg <- paste0(x$has_rejects, " (", round(100 * x$has_rejects / x$n_fits), "%) fits had some steps ",
                  "rejected. Maximum number of rejections was ", x$max_rejects, ".")
    message_list[[2]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[2]] <- data.frame(ok = TRUE, message = "No fits had steps rejected.")
  }

  message_list[[3]] <- data.frame(ok = TRUE, message = paste0("Maximum time was ", x$max_time, " sec."))

  SBC_diagnostic_messages(do.call(rbind, message_list))
}

# For internal use, creates brms backend.
new_SBC_backend_brms <- function(compiled_model,
  args
) {
  require_brms_version("brms backend")

  arg_names_for_stan <- c("chains", "inits", "iter", "warmup", "thin")
  args_for_stan <- args[intersect(names(args), arg_names_for_stan)]
  stan_backend <- sampling_backend_from_stanmodel(compiled_model, args_for_stan)

  structure(list(stan_backend = stan_backend, args = args), class = "SBC_backend_brms")
}

validate_SBC_backend_brms_args <- function(args) {
  if(!is.null(args$algorithm) && args$algorithm != "sampling") {
    stop("Algorithms other than sampling not supported yet")
  }

  unacceptable_params <- c("data", "cores", "empty")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }
}

#' Build a backend based on the `brms` package.
#'
#' @param ... arguments passed to `brm`.
#' @param template_dataset a representative dataset that can be used to generate code.
#' @export
SBC_backend_brms <- function(..., template_dataset) {
  args = list(...)
  validate_SBC_backend_brms_args(args)

  stanmodel <- stanmodel_for_brms(data = template_dataset, ...)

  new_SBC_backend_brms(stanmodel, args)
}

#' Build a brms backend, reusing the compiled model from a previously created `SBC_generator_brms`
#' object.
#'
#' @export
SBC_backend_brms_from_generator <- function(generator, ...) {
  stopifnot(inherits(generator, "SBC_generator_brms"))
  validate_SBC_backend_brms_args(list(...))

  args <- combine_args(generator$args, list(...))
  args$data <- NULL
  args$cores <- NULL
  args$empty <- NULL

  validate_SBC_backend_brms_args(args)


  new_SBC_backend_brms(generator$compiled_model, args)
}

#' @export
SBC_fit.SBC_backend_brms <- function(backend, generated, cores) {
  args_with_data <- backend$args
  args_with_data$data <- generated

  standata <- do.call(brms::make_standata, args_with_data)
  class(standata) <- NULL
  stanfit <- SBC_fit(backend$stan_backend, standata, cores)


  brmsfit <- brmsfit_from_stanfit(stanfit, args_with_data)
  brmsfit
}

#' @export
SBC_fit_to_draws_matrix.brmsfit <- function(fit) {
  posterior::as_draws_matrix(fit$fit)
}

#' @export
SBC_fit_to_diagnostics.brmsfit <- function(fit, fit_output, fit_messages, fit_warnings) {
  SBC_fit_to_diagnostics(fit$fit, fit_output, fit_messages, fit_warnings)
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_brms <- function(backend) {
  object_for_hash <- list(args = backend$args,
                          model_hash =
                            SBC_backend_hash_for_cache(backend$stan_backend))
  rlang::hash(object_for_hash)
}

