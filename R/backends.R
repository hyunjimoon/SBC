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
  if(is.null(args$refresh)) {
    args$refresh <- 0
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
  if(is.null(args$refresh)) {
    args$refresh <- 0
  }
  structure(list(model = model, args = args), class = "cmdstan_sample_SBC_backend")
}

#TODO add SBC_diagnostics generic to extract divergences etc.

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

SBC_fit_to_draws_matrix.brmsfit <- function(fit) {
  posterior::as_draws_matrix(fit$fit)
}

