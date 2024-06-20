#' Backend using bridgesampling to do Bayes factor calculation (and Bayesian Model Averaging) over two candidate
#' models.
#'
#' @param ... passed to [bridgesampling::bridge_sampler()].
#' @export
SBC_backend_bridgesampling <- function(backend_H0, backend_H1, model_var = "model", ...) {
  require_package_version("bridgesampling", version = "1.0", purpose = " to use the bridgesampling SBC backend")
  structure(list(backend_H0 = backend_H0,
                 backend_H1 = backend_H1,
                 model_var = model_var,
                 bridgesampling_args = list(...)),
            class = "SBC_backend_bridgesampling")
}

#' Convert a fit for a given backend into a bridge sampler object.
#' @param ... passed to [bridgesampling::bridge_sampler()].
#' @seealso [bridgesampling::bridge_sampler()]
#' @export
SBC_fit_to_bridge_sampler <- function(backend, fit, generated, ...) {
  UseMethod("SBC_fit_to_bridge_sampler")
}

#' @export
SBC_fit_to_bridge_sampler.default <- function(backend, fit, generated, ...) {
  stop(paste0("To use bridgesampling with backend of class '", class(backend)), "',",
       "you need to implement a corresponding S3 method for SBC_fit_to_bridge_sampler")
}

#' @export
SBC_fit_to_bridge_sampler.SBC_backend_rstan_sample <- function(backend, fit, generated, ...) {
  bridgesampling::bridge_sampler(fit, ...)
}

#' @export
SBC_fit_to_bridge_sampler.SBC_backend_brms <- function(backend, fit, generated, ...) {
  bridgesampling::bridge_sampler(fit, ...)
}


#' @export
SBC_fit.SBC_backend_bridgesampling <- function(backend, generated, cores) {
  fit0 <- SBC_fit(backend$backend_H0, generated, cores)
  fit1 <- SBC_fit(backend$backend_H1, generated, cores)
  bridge_H0 <- do.call("SBC_fit_to_bridge_sampler", c(list(backend$backend_H0, fit0, generated), backend$bridgesampling_args))
  bridge_H1 <- do.call("SBC_fit_to_bridge_sampler", c(list(backend$backend_H1, fit1, generated), backend$bridgesampling_args))
  structure(list(
    fit0 = fit0,
    fit1 = fit1,
    bridge_H0 = bridge_H0,
    bridge_H1 = bridge_H1,
    model_var = backend$model_var
  ), class = "SBC_fit_bridgesampling")
}

SBC_fit_bridgesampling_to_prob1 <- function(fit, log.p = FALSE) {
  log_bf_01 <- bridgesampling::bf(fit$bridge_H0, fit$bridge_H1, log = TRUE)$bf
  prob1 <- plogis(-log_bf_01, log.p = log.p)
  return(prob1)
}

#' @export
SBC_fit_to_draws_matrix.SBC_fit_bridgesampling <- function(fit) {
  draws0 <- posterior::merge_chains(SBC_fit_to_draws_matrix(fit$fit0))
  draws1 <- posterior::merge_chains(SBC_fit_to_draws_matrix(fit$fit1))

  if(posterior::ndraws(draws0) != posterior::ndraws(draws1)) {
    warning("Unequal number of draws for each bridgesampling fit. Will subset to the smaller number.")
    if(posterior::ndraws(draws0) > posterior::ndraws(draws1)) {
      draws0 <- posterior::subset_draws(draws0, draw = 1:posterior::ndraws(draws1))
    } else {
      draws1 <- posterior::subset_draws(draws1, draw = 1:posterior::ndraws(draws0))
    }
  }

  prob1 <- SBC_fit_bridgesampling_to_prob1(fit)

  total_draws <- posterior::ndraws(draws0)

  model_draws <- rbinom(n = total_draws, size = 1, prob = prob1)

  combined_draws <- combine_draws_matrix_for_bf(draws0, draws1, model_draws, model_var = fit$model_var)

  return(combined_draws)
}

#' @export
SBC_posterior_cdf.SBC_fit_bridgesampling <- function(fit, variables) {
  if(fit$model_var %in% names(variables)) {
    prob1 <- SBC_fit_bridgesampling_to_prob1(fit)
    if(variables[fit$model_var] == 0) {
      cdf_low <- 0
      cdf_high <- prob1
    } else if(variables[fit$model_var] == 1) {
      cdf_low <- prob1
      cdf_high <- 1
    } else {
      warning("SBC_posterior_cdf.SBC_fit_bridgesampling expects the model variable to be either 0 or 1")
      cdf_low <- 0
      cdf_high <- 0
    }

    return(data.frame(variable = fit$model_var,
                      cdf_low = cdf_low,
                      cdf_high = cdf_high))
  } else {
    return(NULL)
  }
}

#' @export
SBC_fit_to_diagnostics.SBC_fit_bridgesampling <- function(fit, fit_output, fit_messages, fit_warnings) {
  diags0 <- SBC_fit_to_diagnostics(fit$fit0, fit_output, fit_messages, fit_warnings)
  diags1 <- SBC_fit_to_diagnostics(fit$fit1, fit_output, fit_messages, fit_warnings)

  prob1 <- SBC_fit_bridgesampling_to_prob1(fit)
  log_prob1 <- SBC_fit_bridgesampling_to_prob1(fit, log.p = TRUE)

  percentage_error0 <- bridgesampling::error_measures(fit$bridge_H0)$percentage
  percentage_error1 <- bridgesampling::error_measures(fit$bridge_H1)$percentage
  diags_bs <- data.frame(prob_H1 = prob1, bs_error_H0 = percentage_error0, bs_error_H1 = percentage_error1, log_prob_H1 = log_prob1)

  if(!is.null(diags0)) {
    names(diags0) <- paste0(names(diags0), "_H0")
    diags_bs <- cbind(diags_bs, diags0)
  }

  if(!is.null(diags1)) {
    names(diags1) <- paste0(names(diags1), "_H1")
    diags_bs <- cbind(diags_bs, diags1)
  }

  return(diags_bs)
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_bridgesampling <- function(backend) {
  backend_for_hash <- backend
  backend_for_hash$backend_H0 <- SBC_backend_hash_for_cache(backend$backend_H0)
  backend_for_hash$backend_H1 <- SBC_backend_hash_for_cache(backend$backend_H1)
  rlang::hash(backend_for_hash)
}
