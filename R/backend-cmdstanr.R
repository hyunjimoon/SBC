

#' Backend based on sampling via `cmdstanr`.
#'
#' @param model an object of class `CmdStanModel` (as created by `cmdstanr::cmdstan_model`)
#' @param ... other arguments passed to the `$sample()` method of the model. The `data` and
#'   `parallel_chains` arguments cannot be set this way as they need to be controlled by the SBC
#'   package.
#' @param init_factory an optional function that takes in a dataset and returns a value that
#' can be passed to the `init` argument of `$sample()`. This allows for data-dependent inits.
#' The caching mechanism in [compute_SBC()] ignores the environment of the function, i.e.
#' if the init factory takes values from its environment and those values change between
#' runs, this will not by itself cause cached results to be recomputed.
#' @export
SBC_backend_cmdstan_sample <- function(model, ..., init_factory = NULL) {
  require_cmdstanr_version("cmdstan backend")

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

  if(!is.null(init_factory)) {
    if(any(names(args) == "init")) {
      stop("You cannot specify both init and init_factory")
    }
    stopifnot("init_factory has to be a function" = is.function(init_factory))
    stopifnot("init_factory has to be a function taking a single argument" = length(formalArgs(init_factory)) == 1)
  }
  structure(list(model = model, args = args, init_factory = init_factory), class = "SBC_backend_cmdstan_sample")
}

#' @export
SBC_fit.SBC_backend_cmdstan_sample <- function(backend, generated, cores) {
  our_args <- list(
    data = generated,
    parallel_chains = cores)
  if(!is.null(backend$init_factory)) {
    our_args$init <- backend$init_factory(generated)
  }
  fit <- do.call(backend$model$sample,
                 combine_args(backend$args,
                              our_args))

  if(all(fit$return_codes() != 0)) {
    stop("No chains finished succesfully")
  }

  fit
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_cmdstan_sample <- function(backend) {
  list_for_hash <- list(model = backend$model$code(), args = backend$args)
  # Adding init_factory only when non-null to avoid invalidating old hashes
  if(!is.null(backend$init_factory)) {
    # Explicitly ignoring the environment of the init function,
    # otherwise a refit will be trigerred by a restart of the session
    list_for_hash$init_factory_formals <- formals(backend$init_factory)
    list_for_hash$init_factory_body <- body(backend$init_factory)
  }
  # If there is a user header, add it to the hash
  if(!is.null(backend$model$cpp_options()$USER_HEADER)) {
    header_code <- readLines(backend$model$cpp_options()$USER_HEADER)
    list_for_hash$user_header_code <- header_code
  }
  rlang::hash(list_for_hash)
}


#' @export
SBC_fit_to_draws_matrix.CmdStanMCMC <- function(fit) {
  fit$draws(format = "draws_matrix")
}


#' @export
SBC_fit_to_diagnostics.CmdStanMCMC <- function(fit, fit_output, fit_messages, fit_warnings) {
  diag_summary <- fit$diagnostic_summary()
  res <- data.frame(
    max_chain_time = max(fit$time()$chains[,"total"]),
    n_failed_chains = fit$num_chains() - sum(fit$return_codes() == 0),
    n_divergent = sum(diag_summary$num_divergent),
    n_max_treedepth =  sum(diag_summary$num_max_treedepth),
    min_bfmi = min(diag_summary$ebfmi),
    n_rejects = sum(grepl("reject", fit_messages)) + sum(grepl("reject", fit_warnings))
  )
  res
}

#' @export
SBC_diagnostics_types.CmdStanMCMC <- function(backend) {
  SBC_nuts_diagnostic_types()
}

#' Backend based on variational approximation via `cmdstanr`.
#'
#' @param model an object of class `CmdStanModel` (as created by `cmdstanr::cmdstan_model`)
#' @param n_retries_init number of times to retry the variational fit if the algorithm
#' has trouble initializing (e.g. too many dropped evaluations
#' (see https://discourse.mc-stan.org/t/advi-too-many-dropped-evaluations-even-for-well-behaved-models/24338),
#' or "cannot compute ELBO using the initial variational distribution")
#' @param ... other arguments passed to the `$variational()` method of the model.
#' The `data` argument cannot be set this way as they need to be controlled by the SBC
#'   package.
#' @export
SBC_backend_cmdstan_variational <- function(model, ...,
                                            n_retries_init = 1) {
  require_cmdstanr_version("cmdstan backend")

  stopifnot(inherits(model, "CmdStanModel"))
  if(length(model$exe_file()) == 0) {
    stop("The model has to be already compiled, call $compile() first.")
  }
  args <- list(...)
  unacceptable_params <- c("data")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }

  structure(list(model = model, n_retries_init = n_retries_init, args = args), class = "SBC_backend_cmdstan_variational")
}

stan_variational_elbo_converged <- function(fit_output) {
  any(grepl("ELBO CONVERGED", fit_output))
}

#' @export
SBC_fit.SBC_backend_cmdstan_variational <- function(backend, generated, cores) {
  fit_outputs <- list()
  for(i in 1:backend$n_retries_init) {
    # Need to do my own output capturing as the calling code
    # also captures output and interferes with CmdStanVB$output()
    fit_outputs[[i]] <- capture_all_outputs({
      if(i > 1) {
        cat("==== SBC backend retrying ===== \n")
      }
      fit <- do.call(backend$model$variational,
                     combine_args(backend$args,
                                  list(
                                    data = generated)))
    })

    # Only retry if the error is "too many dropped evaluations" or
    # Cannot compute initial ELBO
    if(fit$return_codes() != 0 &&
       (any(grepl("dropped evaluations.*maximum", fit_outputs[[i]]$messages))
        || any(grepl("Cannot compute ELBO.*initial", fit_outputs[[i]]$messages))
       )
    ) {
      next
    } else {
      break
    }
  }

  reemit_captured(fit_outputs[[i]])

  if(all(fit$return_codes() != 0)) {
    stop("Variational inference did not finish succesfully")
  }

  fit
}


#' @export
SBC_backend_hash_for_cache.SBC_backend_cmdstan_variational <- function(backend) {
  rlang::hash(list(model = backend$model$code(), n_retries_init = backend$n_retries_init, args = backend$args))
}


#' @export
SBC_fit_to_draws_matrix.CmdStanVB <- function(fit) {
  fit$draws(format = "draws_matrix")
}

#' @export
SBC_backend_iid_draws.SBC_backend_cmdstan_variational <- function(backend) {
  TRUE
}


#' @export
SBC_fit_to_diagnostics.CmdStanVB <- function(fit, fit_output, fit_messages, fit_warnings) {
  res <- data.frame(
    elbo_converged = stan_variational_elbo_converged(fit_output),
    n_rejects = sum(grepl("reject", fit_messages)) + sum(grepl("reject", fit_warnings)),
    time = fit$time()$total
  )

  class(res) <- c("SBC_ADVI_diagnostics", class(res))
  res
}

#' @export
SBC_nuts_diagnostic_types.CmdStanVB <- function(backend) {
  SBC_ADVI_diagnostics_types()
}
