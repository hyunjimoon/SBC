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
#' @param fit The fit returned by `SBC_fit`
#' @param fit_output a character string capturing what the backend wrote to stdout
#' @param fit_messages a `data.frame` with columns `type` (type of message as signalled by R,
#'    currently either "message" or "warning") and `message` (the actual text of the message).
#' @return an object that includes diagnostics or NULL, if no diagnostics available.
#' @export
SBC_fit_to_diagnostics <- function(fit, fit_output, fit_messages) {
  UseMethod("SBC_fit_to_diagnostics")
}

#' @export
SBC_fit_to_diagnostics.default <- function(fit, fit_output, fit_messages) {
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

  if(is.null(args$open_progress)) {
    args$open_progress <- FALSE
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
SBC_fit_to_diagnostics.stanfit <- function(fit, fit_output, fit_messages) {
  res <- data.frame(
    max_chain_time = max(rowSums(rstan::get_elapsed_time(fit))),
    n_divergent = rstan::get_num_divergent(fit),
    n_max_treedepth = rstan::get_num_max_treedepth(fit),
    min_bfmi = min(rstan::get_bfmi(fit)),
    n_rejects = sum(grepl("reject", fit_messages$message))
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
    has_treedepth = sum(diagnostics$n_max_treedepth > 0),
    has_rejects = sum(diagnostics$n_rejects > 0),
    max_rejects = max(diagnostics$n_rejects)
  )

  if(!is.null(diagnostics$min_bfmi)) {
    summ$has_low_bfmi = sum(diagnostics$min_bfmi < 0.2)
  }

  if(!is.null(diagnostics$n_failed_chains)) {
    summ$has_failed_chains = sum(diagnostics$n_failed_chains > 0)
  }

  structure(summ, class = "SBC_nuts_diagnostics_summary")
}


get_diagnostics_messages.SBC_nuts_diagnostics <- function(x) {
  get_diagnostics_messages(summary(x))
}

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

get_diagnostics_messages.SBC_nuts_diagnostics_summary <- function(x) {
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
    msg <- paste0(x$has_divergent, " (", round(100 * x$has_divergent / x$n_fits), "%) fits had divergent transitions.")
    message_list[[i]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had divergent transitions.")
  }
  i <- i + 1

  if(x$has_treedepth > 0) {
    msg <- paste0(x$has_treedepth, " (", round(100 * x$has_treedepth / x$n_fits), "%) fits had iterations that saturated max treedepth.")
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
  msg <- get_diagnostics_messages(x)
  print(msg)
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
SBC_fit_to_diagnostics.CmdStanMCMC <- function(fit, fit_output, fit_messages) {
  res <- data.frame(
    max_chain_time = max(fit$time()$chains[,"total"]),
    n_failed_chains = fit$num_chains() - sum(fit$return_codes() == 0),
    n_divergent = sum(fit$sampler_diagnostics()[, , "divergent__"]),
    n_max_treedepth =  sum(fit$sampler_diagnostics()[, , "treedepth__"] == fit$metadata()$max_treedepth),
    n_rejects = sum(grepl("reject", fit_messages$message))
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

#' @export
SBC_fit_to_draws_matrix.brmsfit <- function(fit) {
  posterior::as_draws_matrix(fit$fit)
}

#' @export
SBC_fit_to_diagnostics.brmsfit <- function(fit, fit_output, fit_messages) {
  SBC_fit_to_diagnostics(fit$fit, fit_output, fit_messages)
}
