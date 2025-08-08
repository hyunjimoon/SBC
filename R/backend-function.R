#' A backend that just calls a function.
#'
#' If the function returns a `draws_matrix` object, no other
#' work is necessary to make it work with SBC.
#'
#' @param func the function that will be called in [SBC_fit()]
#' @param generated_arg name of the argument of `func` that will receive the
#' generated data. If `NULL`, data is not passed to the function.
#' @param cores_arg name of the argument of `func` that will receive the number
#' of cores to use. If `NULL`, information on cores is not passed.
#' @param args a (named) list of additional arguments to the function
#' @param iid_draws does the result of the backend have independent identically
#' distribute draws (will be returned by [SBC_backend_iid_draws()] for this backend).
#' @param default_thin_ranks suggested thinning if user does not specify any (will be returned by [SBC_backend_default_thin_ranks()] for this backend).
#' @example examples/backend-function-example.R
#' @export
SBC_backend_function <- function(func, generated_arg = "generated", cores_arg = NULL,
                                 args = list(),
                                 iid_draws = FALSE,
                                 default_thin_ranks = 10) {
  structure(
    list(func = func,
         generated_arg = generated_arg,
         cores_arg = cores_arg,
         iid_draws = iid_draws,
         default_thin_ranks = default_thin_ranks,
         args = args),
    class = "SBC_backend_function"
  )
}

#' @export
SBC_fit.SBC_backend_function <- function(backend, generated, cores) {
  func_args <- list()
  if(!is.null(backend$generated_arg)) {
    func_args[[backend$generated_arg]] <- generated
  }
  if(!is.null(backend$cores_arg)) {
    func_args[[backend$cores_arg]] <- cores
  }

  func_args <- c(func_args, backend$args)
  do.call(backend$func, func_args)
}

#' @export
SBC_backend_iid_draws.SBC_backend_function <- function(backend) {
  backend$iid_draws
}

#' @export
SBC_backend_default_thin_ranks.SBC_backend_function <- function(backend) {
  backend$default_thin_ranks
}
