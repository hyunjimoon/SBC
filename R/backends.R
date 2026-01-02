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

#' Use a fitted model to obtain a posterior CDF at a given point.
#'
#' Backends that represent explicit posterior distributions (i.e. not just samples)
#' can implement this S3 generic to evaluate the CDF of that distribution at
#' the simulated values.
#'
#' A backend may choose to return explicit CDF only for some model parameters.
#' return from [SBC_fit()].
#'
#' @param fit an object returned by the [SBC_fit()] method.
#' @param variables a named vector of values at which to evaluate the CDF
#' @return either `NULL` (the default implementation) or a `data.frame`
#' with a row for each variable. The columns must include:
#' `variable` (variable name), and either `cdf` (if there are no ties) or the
#' pair `cdf_low` and `cdf_high` indicate the possible tied values for the CDF.
#'
#' @export
SBC_posterior_cdf <- function(fit, variables) {
  UseMethod("SBC_posterior_cdf")
}

#' @rdname SBC_posterior_cdf
#' @export
SBC_posterior_cdf.default <- function(fit, variables) {
  NULL
}

#' @return A [derived_quantities()] object or NULL if there are none
#' @export
SBC_fit_specific_dquants <- function(fit) {
  UseMethod("SBC_fit_specific_dquants")
}


#' @rdname SBC_fit_specific_dquants
#' @export
SBC_fit_specific_dquants.default <- function(fit) {
  NULL
}

#' Validate the results of a SBC_posterior_cdf call.
#' @keywords internal
validate_cdf_df <- function(cdf_df, valid_variables) {
  if(!is.null(cdf_df)) {
    if(!is.data.frame(cdf_df)) {
      stop("Result of SBC_posterior_cdf must either be NULL or a data.frame")
    }
    if(!("variable" %in% names(cdf_df))) {
      stop("Result of SBC_posterior_cdf must have a 'variable' column.")
    }
    if(!(("cdf_low" %in% names(cdf_df)) || !("cdf_high" %in% names(cdf_df)))) {
      if(!("cdf" %in% names(cdf_df))) {
        stop("Result of SBC_posterior_cdf must either have a 'cdf' column or both a 'cdf_low' and a 'cdf_high' column.")
      } else {
        cdf_df$cdf_low <- cdf_df$cdf
        cdf_df$cdf_high <- cdf_df$cdf
      }
    }
    unknown_vars <- cdf_df$variable[!(cdf_df$variable %in% valid_variables)]
    if(length(unknown_vars) > 0) {
      stop(paste0("Result of SBC_posterior_cdf refers to unknown variables: ",
                  paste0(unknown_vars, ", ")))
    }
  }
  return(cdf_df)
}

#' S3 generic to get backend-specific diagnostics.
#'
#' The diagnostics object has to be a `data.frame` but may
#' inherit additional classes - in particular it may be useful
#' for the returning object to implement [diagnostic_types()].
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

#' @export
diagnostic_types <- function(diags) {
  UseMethod("diagnostic_types")
}

#' @export
diagnostic_types.default <- function(diags) {
  list()
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
