#' @export
SBC_fit <- function(backend, generated, cores) {
  UseMethod("SBC_fit")
}

#' @export
cmdstan_sample_SBC_backend <- function(model, ...) {
  stopifnot(inherits(model, "CmdStanModel"))
  #TODO: check data, cores not in args
  structure(list(model = model, args = list(...)), class = "cmdstan_sample_SBC_backend")
}

draws_rvars_to_standata <- function(x) {
  stopifnot(posterior::ndraws(x) == 1)
  lapply(x, FUN = function(x_rvar) {
    res <- draws_of(x_rvar, with_chains = FALSE)
    #TODO figure out how to distinguish between scalar and array of size 1
    if(identical(dim(x_rvar), 1L)) {
      as.numeric(res[1,])
    } else {
      dim(res) <- dim(res)[2:length(dim(res))]
      res
    }
  })
}

#TODO add SBC_diagnostics generic to extract divergences etc.

SBC_fit.cmdstan_sample_SBC_backend <- function(backend, generated, cores) {
  do.call(backend$model$sample,
          c(list(refresh = 0),
            backend$args,
            list(
              data = generated,
              parallel_chains = cores)))
}

#' @export
SBC_fit_to_draws_rvars <- function(fit) {
  UseMethod("SBC_fit_to_draws_rvars")
}

SBC_fit_to_draws_rvars.CmdStanMCMC <- function(fit) {
  posterior::as_draws_rvars(fit$draws())
}
