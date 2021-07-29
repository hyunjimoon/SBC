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
  x
}

SBC_fit.cmdstan_sample_SBC_backend <- function(backend, generated, cores) {
  backend$model$sample(c(backend$args),
                   list(
                     data = draws_rvars_to_standata(generated),
                     cores = cores))
}

#' @export
SBC_fit_to_draws_rvars <- function(fit) {
  UseMethod("SBC_fit_to_draws_rvars")
}

SBC_fit_to_draws_rvars.CmdStanMCMC(fit) {
  fit$draws(format = "draws_rvars")
}

