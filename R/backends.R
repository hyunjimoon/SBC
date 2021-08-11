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

#' @export
SBC_fit_to_diagnostics <- function(fit) {
  UseMethod("SBC_fit_to_diagnostics")
}

#' @export
SBC_fit_to_diagnostics.default <- function(fit) {
  NULL
}

#' @export
rstan_sample_SBC_backend <- function(model, ...) {
  stopifnot(inherits(model, "stanmodel"))
  args <- list(...)
  unacceptable_params <- c("data", "cores")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }
  if(is.null(args$refresh)) {
    args$refresh <- 0
  }
  structure(list(model = model, args = args), class = "rstan_sample_SBC_backend")
}


#' @export
SBC_fit.rstan_sample_SBC_backend <- function(backend, generated, cores) {
  do.call(rstan::sampling,
          combine_args(list(object = backend$model,
                 data = generated,
                 cores = cores),
            backend$args
            ))
}

#' @export
SBC_fit_to_diagnostics.stanfit <- function(fit) {
  res <- data.frame(
    max_chain_time = max(rowSums(rstan::get_elapsed_time(fit))),
    n_divergent = rstan::get_num_divergent(fit),
    n_max_treedepth = rstan::get_num_max_treedepth(fit),
    min_bfmi = min(rstan::get_bfmi(fit))
  )

  class(res) <- c("SBC_nuts_diagnostics", class(res))
  res
}

#' @export
summary.SBC_nuts_diagnostics <- function(diagnostics) {
  summ <- list(
    n_fits = nrow(diagnostics),
    max_chain_time = max(diagnostics$max_chain_time),
    has_divergent = sum(diagnostics$n_divergent > 0),
    has_treedepth = sum(diagnostics$n_max_treedepth > 0)
  )

  if(!is.null(diagnostics$min_bfmi)) {
    summ$has_low_bfmi = sum(diagnostics$min_bfmi < 0.2)
  }

  if(!is.null(diagnostics$n_failed_chains)) {
    summ$has_failed_chains = sum(diagnostics$n_failed_chains > 0)
  }

  structure(summ, class = "SBC_nuts_diagnostics_summary")
}

#' @export
print.SBC_nuts_diagnostics_summary <- function(x) {
  if(!is.null(x$has_failed_chains)) {
    if(x$has_failed_chains > 0) {
      cat(" - ", x$has_failed_chains, " (", round(100 * x$has_failed_chains / x$n_fits), "%) fits had some failed chains.\n", sep = "")
    } else {
      cat(" - No fits had failed chains.\n")
    }
  }

  if(x$has_divergent > 0) {
    cat(" - ", x$has_divergent, " (", round(100 * x$has_divergent / x$n_fits), "%) fits had divergent transitions.\n", sep = "")
  } else {
    cat(" - No fits had divergent transitions.\n")
  }
  if(x$has_treedepth > 0) {
    cat(" - ", x$has_treedepth, " (", round(100 * x$has_treedepth / x$n_fits), "%) fits had iterations that saturated max treedepth.\n", sep = "")
  } else {
    cat(" - No fits had iterations that saturated max treedepth.\n")
  }

  if(!is.null(x$has_low_bfmi)) {
    if(x$has_low_bfmi > 0) {
      cat(" - ", x$has_low_bfmi, " (", round(100 * x$has_low_bfmi / x$n_fits), "%) fits had low BFMI.\n", sep = "")
    } else {
      cat(" - No fits had low BFMI.\n")
    }
  }
  invisible(x)
}

#' @export
cmdstan_sample_SBC_backend <- function(model, ...) {
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
  if(is.null(args$refresh)) {
    args$refresh <- 0
  }
  structure(list(model = model, args = args), class = "cmdstan_sample_SBC_backend")
}

#' @export
SBC_fit.cmdstan_sample_SBC_backend <- function(backend, generated, cores) {
  do.call(backend$model$sample,
          combine_args(backend$args,
            list(
              data = generated,
              parallel_chains = cores)))
}



#' @export
SBC_fit_to_draws_matrix.CmdStanMCMC <- function(fit) {
  fit$draws(format = "draws_matrix")
}


#' @export
SBC_fit_to_diagnostics.CmdStanMCMC <- function(fit) {
  res <- data.frame(
    max_chain_time = max(fit$time()$chains[,"total"]),
    n_failed_chains = fit$num_chains() - sum(fit$return_codes() == 0),
    n_divergent = sum(fit$sampler_diagnostics()[, , "divergent__"]),
    n_max_treedepth =  sum(fit$sampler_diagnostics()[, , "treedepth__"] == fit$metadata()$max_treedepth)
    #
  ) # TODO: add min_bfmi once https://github.com/stan-dev/cmdstanr/pull/500/ is merged
  class(res) <- c("SBC_nuts_diagnostics", class(res))
  res
}


new_brms_SBC_backend <- function(compiled_model,
  args
) {

  arg_names_for_stan <- c("chains", "inits", "iter", "warmup", "thin")
  args_for_stan <- args[intersect(names(args), arg_names_for_stan)]
  stan_backend <- sampling_backend_from_stanmodel(compiled_model, args_for_stan)

  structure(list(stan_backend = stan_backend, args = args), class = "brms_SBC_backend")
}

validate_brms_SBC_backend_args <- function(args) {
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

#' Build a brms backend.
#'
#' @param ... arguments passed to `brm`.
#' @param template_dataset a representative dataset that can be used to generate code.
#' @export
brms_SBC_backend <- function(..., template_dataset) {
  args = list(...)
  validate_brms_SBC_backend_args(args)

  stanmodel <- stanmodel_for_brms(data = template_dataset, ...)

  new_brms_SBC_backend(stanmodel, args)
}

#' Build a brms backend, reusing the compiled model from a previously created generator.
#' @export
brms_SBC_backend_from_generator <- function(generator, ...) {
  stopifnot(inherits(generator, "brms_SBC_generator"))
  validate_brms_SBC_backend_args(list(...))

  args <- combine_args(generator$args, list(...))

  if(!is.null(args$algorithm) && args$algorithm != "sampling") {
    stop("Algorithms other than sampling not supported yet")
  }

  new_brms_SBC_backend(generator$compiled_model, args)
}

#' @export
SBC_fit.brms_SBC_backend <- function(backend, generated, cores) {
  args_with_data <- backend$args
  args_with_data$data <- generated

  standata <- do.call(brms::make_standata, args_with_data)
  class(standata) <- NULL
  stanfit <- SBC_fit(backend$stan_backend, standata, cores)


  brmsfit <- brmsfit_from_stanfit(stanfit, args_with_data)
  brmsfit
}

SBC_fit_to_draws_matrix.brmsfit <- function(fit) {
  posterior::as_draws_matrix(fit$fit)
}

