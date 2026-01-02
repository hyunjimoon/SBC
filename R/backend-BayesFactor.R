#' A backend using [BayesFactor::lmBF()]
#'
#' @param ... all parameters are passed directly to [BayesFactor::lmBF()].
#' @export
SBC_backend_lmBF <- function(..., iterations = 2000) {
  SBC:::require_package_version("BayesFactor", version = "0.9", purpose = " to use the lmBF SBC backend")
  unacceptable_params <- c("data", "posterior")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }
  args_list <- list(...)
  args_list$iterations <- iterations
  structure(list(
    args = args_list
  ), class = "SBC_backend_lmBF")
}

#' @export
SBC_fit.SBC_backend_lmBF <- function(backend, generated, cores) {
  args <- backend$args
  args$data <- generated
  args$posterior <- TRUE
  posterior <- do.call(BayesFactor::lmBF, args)

  args$posterior <- FALSE
  bf <- do.call(BayesFactor::lmBF, args)

  structure(list(posterior = posterior, bf = bf), class = "SBC_fit_lmBF")
}

#' @export
SBC_fit_to_draws_matrix.SBC_fit_lmBF <- function(fit) {
  posterior::as_draws_matrix(fit$posterior)
}

#' A backend for Bayes factor workflow with [BayesFactor::extractBF].
#'
#'
#' @export
SBC_backend_extractBF_comparison <- function(backend_H0, backend_H1, model_var = "model") {
  SBC:::require_package_version("BayesFactor", version = "0.9", purpose = " to use the extractBF_comparison SBC backend")
  structure(list(backend_H0 = backend_H0,
                 backend_H1 = backend_H1,
                 model_var = model_var),
            class = "SBC_backend_extractBF_comparison")
}

#' A fits need to implement this S3 generic to be useful for the
#' [SBC_backend_extractBF_comparison()] backend.
#' @export
SBC_fit_to_BFBayesFactor <- function(fit) {
  UseMethod("SBC_fit_to_BFBayesFactor", fit)
}

#' @export
SBC_fit_to_BFBayesFactor.SBC_fit_lmBF <- function(fit) {
  fit$bf
}



#' @export
SBC_fit.SBC_backend_extractBF_comparison <- function(backend, generated, cores) {
  fit0 <- SBC_fit(backend$backend_H0, generated, cores)
  fit1 <- SBC_fit(backend$backend_H1, generated, cores)
  logbf10 <- BayesFactor::extractBF(SBC_fit_to_BFBayesFactor(fit1)/SBC_fit_to_BFBayesFactor(fit0), logbf = TRUE)
  structure(list(
    fit0 = fit0,
    fit1 = fit1,
    logbf10 = logbf10,
    model_var = backend$model_var
  ), class = "SBC_fit_extractBF_comparison")
}

SBC_fit_extractBF_comparison_to_prob1 <- function(fit, log.p = FALSE) {
  logbf10 <- fit$logbf10$bf
  prob1 <- plogis(logbf10, log.p = log.p)
  return(prob1)
}

#' @export
SBC_fit_to_draws_matrix.SBC_fit_extractBF_comparison <- function(fit) {
  draws0 <- posterior::merge_chains(SBC_fit_to_draws_matrix(fit$fit0))
  draws1 <- posterior::merge_chains(SBC_fit_to_draws_matrix(fit$fit1))

  if(posterior::ndraws(draws0) != posterior::ndraws(draws1)) {
    warning("Unequal number of draws for each extractBF fit. Will subset to the smaller number.")
    if(posterior::ndraws(draws0) > posterior::ndraws(draws1)) {
      draws0 <- posterior::subset_draws(draws0, draw = 1:posterior::ndraws(draws1))
    } else {
      draws1 <- posterior::subset_draws(draws1, draw = 1:posterior::ndraws(draws0))
    }
  }

  prob1 <- SBC_fit_extractBF_comparison_to_prob1(fit)

  total_draws <- posterior::ndraws(draws0)

  model_draws <- rbinom(n = total_draws, size = 1, prob = prob1)

  combined_draws <- SBC:::combine_draws_matrix_for_bf(list(draws0, draws1), model_draws, model_var = fit$model_var)

  return(combined_draws)
}

#' @export
SBC_posterior_cdf.SBC_fit_extractBF_comparison <- function(fit, variables) {
  if(fit$model_var %in% names(variables)) {
    prob1 <- SBC_fit_extractBF_comparison_to_prob1(fit)
    return(binary_to_cdf(fit$model_var, prob1, variables[fit$model_var]))
  } else {
    return(NULL)
  }
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_extractBF_comparison <- function(backend) {
  backend_for_hash <- backend
  backend_for_hash$backend_H0 <- SBC_backend_hash_for_cache(backend$backend_H0)
  backend_for_hash$backend_H1 <- SBC_backend_hash_for_cache(backend$backend_H1)
  rlang::hash(backend_for_hash)
}
