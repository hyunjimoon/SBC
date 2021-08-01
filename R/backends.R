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
SBC_fit_to_draws_matrix <- function(fit) {
  UseMethod("SBC_fit_to_draws_matrix")
}

SBC_fit_to_draws_matrix.default <- function(fit) {
  posterior::as_draws_matrix(fit)
}

SBC_fit_to_draws_matrix.CmdStanMCMC <- function(fit) {
  fit$draws(format = "draws_matrix")
}


SBC_fit_to_draws_matrix.brmsfit <- function(fit) {
  posterior::as_draws_matrix(fit$fit)
}

new_brms_SBC_backend <- function(compiled_model,
  args
) {

  arg_names_for_stan <- c("chains", "inits", "iter", "warmup", "thin")
  args_for_stan <- args[intersect(names(args), arg_names_for_stan)]
  stan_backend <- sampling_backend_from_stanmodel(compiled_model, args_for_stan)

  structure(list(stan_backend = stan_backend, args = args), class = "brms_SBC_backend")
}

#' Build a brms backend.
#'
#' @param ... arguments passed to `brm`.
#' @export
brms_SBC_backend <- function(...) {
  args = list(...)
  if(!is.null(args$algorithm) && args$algorithm != "sampling") {
    stop("Algorithms other than sampling not supported yet")
  }

  stanmodel <- stanmodel_for_brms(...)

  new_brms_SBC_backend(stanmodel, args)
}

#' Build a brms backend, reusing the compiled model from a previously created generator.
#' @export
brms_SBC_backend_from_generator <- function(generator, ...) {
  stopifnot(inherits(generator, "brms_SBC_generator"))
  args <- c(generator$args, list(...))

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
