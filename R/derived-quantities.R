#' Create a definition of derived quantities evaluated in R.
#'
#' When the expression contains non-library functions/objects, and parallel processing
#' is enabled, those must be
#' named in the `.globals` parameter (hopefully we'll be able to detect those
#' automatically in the future). Note that [recompute_SBC_statistics()] currently
#' does not use parallel processing, so `.globals` don't need to be set.
#'
#' @param ... named expressions representing the quantitites
#' @param .globals A list of names of objects that are defined
#' in the global environment and need to present for the gen. quants. to evaluate.
#' It is added to the `globals` argument to [future::future()], to make those
#' objects available on all workers.
#' @param .var_attributes a [var_attributes()] object providing attributes for
#' the derived quantities, if necessary
#' @examples
#'# Derived quantity computing the total log likelihood of a normal distribution
#'# with known sd = 1
#'normal_lpdf <- function(y, mu, sigma) {
#'  sum(dnorm(y, mean = mu, sd = sigma, log = TRUE))
#'}
#'
#'# Note the use of .globals to make the normal_lpdf function available
#'# within the expression
#'log_lik_dq <- derived_quantities(log_lik = normal_lpdf(y, mu, 1),
#'                                 .globals = "normal_lpdf" )
#'
#' @export
derived_quantities <- function(..., .var_attributes = NULL, .globals = list()) {
  structure(rlang::enquos(..., .named = TRUE),
            class = "SBC_derived_quantities",
            globals = .globals,
            var_attributes = .var_attributes
            )
}

#' @title Validate a definition of derived quantities evaluated in R.
#' @export
validate_derived_quantities <- function(x) {
  # Backwards compatibility
  if(inherits(x, "SBC_generated_quantities")) {
    class(x) <- "SBC_derived_quantities"
  }
  stopifnot(inherits(x, "SBC_derived_quantities"))
  invisible(x)
}

#' @title Combine two lists of derived quantities
#' @export
bind_derived_quantities <- function(dq1, dq2) {
  if(is.null(dq1)) {
    return(dq2)
  } else if(is.null(dq2)) {
    return(dq1)
  }
  validate_derived_quantities(dq1)
  validate_derived_quantities(dq2)
  structure(c(dq1, dq2),
            class = "SBC_derived_quantities",
            globals = bind_globals(attr(dq1, "globals"), attr(dq2, "globals")),
            var_attributes = combine_var_attributes(attr(dq1, "var_attributes"), attr(dq2, "var_attributes")))
}

#'@title Compute derived quantities based on given data and posterior draws.
#'@param gen_quants Deprecated, use `dquants`
#'@export
compute_dquants <- function(draws, generated, dquants, gen_quants = NULL) {
  if(!is.null(gen_quants)) {
    warning("gen_quants argument is deprecated, use dquants")
    if(rlang::is_missing(dquants)) {
      dquants <- gen_quants
    }
  }
  dquants <- validate_derived_quantities(dquants)
  draws_rv <- posterior::as_draws_rvars(draws)

  draws_env <- list2env(draws_rv)
  if(!is.null(generated)) {
    if(!is.list(generated)) {
      stop("compute_dquants assumes that generated is a list, but this is not the case")
    }
    generated_env <- list2env(generated, parent = draws_env)

    data_mask <- rlang::new_data_mask(bottom = generated_env, top = draws_env)
  } else {
    data_mask <- rlang::new_data_mask(bottom = draws_env)
  }

  eval_func <- function(dq) {
    # Wrap the expression in `rdo` which will mostly do what we need
    # all the tricks are just to have the correct environment when we need it
    wrapped_dq <- rlang::new_quosure(rlang::expr(posterior::rdo(!!rlang::get_expr(dq))), rlang::get_env(dq))
    rlang::eval_tidy(wrapped_dq, data = data_mask)
  }
  rvars <- lapply(dquants, FUN = eval_func)
  do.call(posterior::draws_rvars, rvars)
}

#' Get the [var_attributes()] object associated with the derived quantities
#' @export
dquants_var_attributes <- function(dquants) {
  return(attr(dquants, "var_attributes", exact = TRUE))
}

#' @title Create a definition of derived quantities evaluated in R.
#' @description Delegates directly to `derived_quantities()`.
#'
#' @name generated_quantities-deprecated
#' @seealso \code{\link{SBC-deprecated}}
#' @keywords internal
NULL

#' @rdname SBC-deprecated
#' @section \code{generated_quantities}:
#' Instead of \code{generated_quantities}, use \code{\link{derived_quantities}}.
#'
#' @export
generated_quantities <- function(...) {
  warning("generated_quantities() is deprecated, use derived_quantities instead.")
  derived_quantities(...)
}

#' @title Validate a definition of derived quantities evaluated in R.
#' @description Delegates directly to `validate_derived_quantities()`.
#'
#' @name generated_quantities-deprecated
#' @seealso \code{\link{SBC-deprecated}}
#' @keywords internal
NULL

#' @rdname SBC-deprecated
#' @section \code{validate_generated_quantities}:
#' Instead of \code{validate_generated_quantities}, use \code{\link{validate_derived_quantities}}.
#'
#' @export
validate_generated_quantities <- function(...) {
  warning("generated_quantities() is deprecated, use validate_derived_quantities instead.")
  validate_derived_quantities(...)
}

#' @title Combine two lists of derived quantities
#' @description Delegates directly to `bind_derived_quantities()`.
#'
#' @name bind_generated_quantities-deprecated
#' @seealso \code{\link{SBC-deprecated}}
#' @keywords internal
NULL

#' @rdname SBC-deprecated
#' @section \code{bind_generated_quantities}:
#' Instead of \code{bind_generated_quantities}, use \code{\link{bind_derived_quantities}}.
#'
#' @export
bind_generated_quantities <- function(...) {
  warning("bind_generated_quantities() is deprecated, use bind_derived_quantities instead.")
  bind_derived_quantities(...)
}

#'@title Compute derived quantities based on given data and posterior draws.
#' @description Delegates directly to `compute_dquants()`.
#'
#' @name compute_gen_quants-deprecated
#' @seealso \code{\link{SBC-deprecated}}
#' @keywords internal
NULL

#' @rdname SBC-deprecated
#' @section \code{compute_gen_quants}:
#' Instead of \code{compute_gen_quants}, use \code{\link{compute_dquants}}.
#'
#' @export
compute_gen_quants <- function(...) {
  warning("compute_gen_quants() is deprecated, use compute_dquants() instead.")
  compute_dquants(...)
}
