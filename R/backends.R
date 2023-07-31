#' S3 generic using backend to fit a model to data.
#'
#' Needs to be implemented by all backends.
#' All implementations have to return an object for which you can safely
#' call [SBC_fit_to_draws_matrix()] and get some draws.
#' If that's not possible an error should be raised.
#' @export
SBC_fit <- function(backend, generated, cores) {
  UseMethod("SBC_fit")
}

#' S3 generic converting a fitted model to a `draws_matrix` object.
#'
#' Needs to be implemented for all types of objects the backend can
#' return from [SBC_fit()]. Default implementation just calls,
#' [posterior::as_draws_matrix()], so if the fit object already supports
#' this, it will work out of the box.
#' @export
SBC_fit_to_draws_matrix <- function(fit) {
  UseMethod("SBC_fit_to_draws_matrix")
}

#' @rdname SBC_fit_to_draws_matrix
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

#' S3 generic to let backends signal that they produced independent draws.
#'
#' Most backends (e.g. those based on variatns of MCMC) don't produce
#' independent draws and thus diagnostics like Rhat and ESS are important
#' and draws may need thinning. Backends that already produce independent
#' draws (e.g. ADVI/optimizing) can implement this method to return `TRUE`
#' to signal this is the case. If this method returns `TRUE`, ESS and Rhat will
#' always attain their best possible values and [SBC_backend_default_thin_ranks()]
#' will return `1`.
#'  The default implementation returns `FALSE`.
#' @param backend to check
#' @export
SBC_backend_iid_draws <- function(backend) {
  UseMethod("SBC_backend_iid_draws")
}

#' @rdname SBC_backend_iid_draws
#' @export
SBC_backend_iid_draws.default <- function(backend) {
  FALSE
}

#' S3 generic to get backend-specific default thinning for rank computation.
#'
#' The default implementation plays it relatively safe and returns 10, unless
#' [SBC_backend_iid_draws()] returns `TRUE` in which case it returns 1.
#'
#' @export
SBC_backend_default_thin_ranks <- function(backend) {
  UseMethod("SBC_backend_default_thin_ranks")
}

#' @rdname SBC_backend_default_thin_ranks
#' @export
SBC_backend_default_thin_ranks.default <- function(backend) {
  if(SBC_backend_iid_draws(backend)) {
    1
  } else {
    10
  }
}
