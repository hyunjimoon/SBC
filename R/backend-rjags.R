#' Create a JAGS backend using `rjags`
#'
#' @param file model file or connection to model code (passed to [rjags::jags.model()])
#' @param n.iter number of iterations for sampling (passed to [rjags::coda.samples())
#' @param n.burnin number of iterations used for burnin
#' @param variable.names names of variables to monitor (passed to [rjags::coda.samples()])
#' @param thin thinning (passed to [rjags::coda.samples()])
#' @param na.rm whether to omit variables containing NA (passed to [rjags::coda.samples()])
#' @param ... additional arguments passed to [rjags::jags.model()]
#' @export
SBC_backend_rjags <- function(file, n.iter, n.burnin, variable.names, thin = 1, na.rm = TRUE, ...) {
  args = list(...)
  if(any(names(args) == "data")) {
    stop(paste0("Argument 'data' cannot be provided when defining a backend",
                " as it needs to be set by the SBC package"))
  }

  structure(list(file = file,
                 n.iter = n.iter,
                 variable.names = variable.names,
                 n.burnin = n.burnin,
                 thin = thin,
                 na.rm = na.rm,
                 args = args), class = "SBC_backend_rjags")
}


#' @export
SBC_fit.SBC_backend_rjags <- function(backend, generated, cores) {
  args_all <- c(list(file = backend$file, data = generated), backend$args)

  model <- do.call(rjags::jags.model, args_all)
  if(backend$n.burnin > 0) {
    stats::update(model, n.iter = backend$n.burnin, progress.bar = "none")
  }
  samples <- rjags::coda.samples(model,
                          variable.names = backend$variable.names,
                          n.iter = backend$n.iter,
                          thin = backend$thin,
                          na.rm = backend$na.rm,
                          progress.bar = "none")

  structure(list(model = model, samples = samples),
            class = "SBC_rjags_fit")
}


#' @export
SBC_fit_to_draws_matrix.SBC_rjags_fit <- function(fit) {
  posterior::as_draws_matrix(fit$samples)
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_rjags <- function(backend) {
  model_code <- readLines(backend$file)

  backend_for_cache <- backend
  backend_for_cache$file <- NULL
  backend_for_cache$model_code <- model_code
  rlang::hash(backend_for_cache)
}
