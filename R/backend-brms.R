

# For internal use, creates brms backend.
new_SBC_backend_brms <- function(compiled_model,
                                 args
) {
  require_brms_version("brms backend")

  if(args$algorithm == "sampling") {
    arg_names_for_stan <- c("chains", "inits", "init", "iter", "warmup", "thin")
  } else if(args$algorithm == "meanfield") {
    arg_names_for_stan <- c("inits", "init") # possibly more valid args
  }
  args_for_stan <- args[intersect(names(args), arg_names_for_stan)]

  args_for_stan_renames <- c("inits" = "init")
  for(i in 1:length(args_for_stan_renames)) {
    orig <- names(args_for_stan_renames)[i]
    new <- args_for_stan_renames[i]
    if(!is.null(args_for_stan[[orig]])) {
      args_for_stan[[new]] <- args_for_stan[[orig]]
      args_for_stan[[orig]] <- NULL
    }
  }
  if(args$algorithm == "sampling") {
    stan_backend <- sampling_backend_from_stanmodel(compiled_model, args_for_stan)
  } else if(args$algorithm == "meanfield") {
    stan_backend <- variational_backend_from_stanmodel(compiled_model, args_for_stan)
  }

  structure(list(stan_backend = stan_backend, args = args), class = "SBC_backend_brms")
}

validate_SBC_backend_brms_args <- function(args) {
  if(!is.null(args$algorithm) && args$algorithm != "sampling" && args$algorithm != "meanfield") {
    stop("Algorithms other than sampling and meanfield not supported yet. Comment on https://github.com/hyunjimoon/SBC/issues/91 to express your interest.")
  }

  unacceptable_params <- c("data", "cores", "empty")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }
}

#' Build a backend based on the `brms` package.
#'
#' @param ... arguments passed to `brm`.
#' @param template_data a representative value for the `data` argument in `brm`
#'    that can be used to generate code.
#' @param template_dataset DEPRECATED. Use `template_data`
#' @param out_stan_file A filename for the generated Stan code. Useful for
#'    debugging and for avoiding unnecessary recompilations.
#' @export
SBC_backend_brms <- function(..., template_data, out_stan_file = NULL, template_dataset = NULL) {
  if(!is.null(template_dataset)) {
    warning("Argument 'template_dataset' is deprecated, use 'template_data' instead")
    if(missing(template_data)) {
      template_data <- template_dataset
    }
  }
  args = list(...)
  validate_SBC_backend_brms_args(args)

  stanmodel <- stanmodel_for_brms(data = template_data, out_stan_file = out_stan_file, ...)

  new_SBC_backend_brms(stanmodel, args)
}

#' Build a brms backend, reusing the compiled model from a previously created `SBC_generator_brms`
#' object.
#'
#' @export
SBC_backend_brms_from_generator <- function(generator, ...) {
  stopifnot(inherits(generator, "SBC_generator_brms"))
  validate_SBC_backend_brms_args(list(...))

  args <- combine_args(generator$args, list(...))
  args$data <- NULL
  args$cores <- NULL
  args$empty <- NULL

  validate_SBC_backend_brms_args(args)


  new_SBC_backend_brms(generator$compiled_model, args)
}

#' @export
SBC_fit.SBC_backend_brms <- function(backend, generated, cores) {
  args_with_data <- backend$args
  args_with_data$data <- generated

  standata <- do.call(brms::make_standata, args_with_data)
  class(standata) <- NULL
  stanfit <- SBC_fit(backend$stan_backend, standata, cores)


  brmsfit <- brmsfit_from_stanfit(stanfit, args_with_data)
  brmsfit
}

#' @export
SBC_fit_to_draws_matrix.brmsfit <- function(fit) {
  posterior::as_draws_matrix(fit$fit)
}

#' @export
SBC_fit_to_diagnostics.brmsfit <- function(fit, fit_output, fit_messages, fit_warnings) {
  SBC_fit_to_diagnostics(fit$fit, fit_output, fit_messages, fit_warnings)
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_brms <- function(backend) {
  object_for_hash <- list(args = backend$args,
                          model_hash =
                            SBC_backend_hash_for_cache(backend$stan_backend))
  rlang::hash(object_for_hash)
}

