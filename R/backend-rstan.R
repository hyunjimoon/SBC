#' SBC backend using the `sampling` method from `rstan`.
#'
#' @param model a `stanmodel` object (created via `rstan::stan_model`)
#' @param ... other arguments passed to `sampling` (number of iterations, ...).
#'   Arguments `data` and `cores` cannot be set this way as they need to be
#'   controlled by the package.
#' @export
SBC_backend_rstan_sample <- function(model, ...) {
  stopifnot(inherits(model, "stanmodel"))
  args <- list(...)
  unacceptable_params <- c("data", "cores")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }

  if(is.null(args$open_progress)) {
    args$open_progress <- FALSE
  }
  structure(list(model = model, args = args), class = "SBC_backend_rstan_sample")
}


#' @export
SBC_fit.SBC_backend_rstan_sample <- function(backend, generated, cores) {
  fit <- do.call(rstan::sampling,
                 combine_args(list(object = backend$model,
                                   data = generated,
                                   ## TODO: Forcing a single core until we can capture output with multiple cores
                                   ## https://discourse.mc-stan.org/t/capturing-warnings-rejects-from-rstan-with-multiple-cores/23976
                                   cores = 1),
                              backend$args
                 ))

  if(fit@mode != 0) {
    stop("Fit does not contain draws.")
  }

  fit
}

#' @export
SBC_fit_to_diagnostics.stanfit <- function(fit, fit_output, fit_messages, fit_warnings) {
  res <- data.frame(
    max_chain_time = max(rowSums(rstan::get_elapsed_time(fit))),
    n_divergent = rstan::get_num_divergent(fit),
    n_max_treedepth = rstan::get_num_max_treedepth(fit),
    min_bfmi = min(rstan::get_bfmi(fit)),
    n_rejects = sum(grepl("reject", fit_messages)) + sum(grepl("reject", fit_warnings))
  )

  class(res) <- c("SBC_nuts_diagnostics", class(res))
  res
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_rstan_sample <- function(backend) {
  rlang::hash(list(model = backend$model@model_code, args = backend$args))
}


#' SBC backend using the `optimizing` method from `rstan`.
#'
#' @param model a `stanmodel` object (created via `rstan::stan_model`)
#' @param ... other arguments passed to `optimizing` (number of iterations, ...).
#'   Argument `data` cannot be set this way as they need to be
#'   controlled by the package.
#' @param n_retries_hessian the number of times the backend is allow to retry optimization
#' (with different seeed) to produce a usable Hessian that can produce draws. In some cases,
#' the Hessian may be numerically unstable and not be positive definite.
#' @export
SBC_backend_rstan_optimizing <- function(model, ..., n_retries_hessian = 1) {
  stopifnot(inherits(model, "stanmodel"))
  n_retries_hessian <- as.integer(n_retries_hessian)
  stopifnot(length(n_retries_hessian) == 1)
  stopifnot(n_retries_hessian > 0)

  args <- list(...)
  unacceptable_params <- c("data", "hessian")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }

  args$hessian <- TRUE
  if(is.null(args$draws)) {
    args$draws <- 1000
  } else if(args$draws <= 1) {
    stop("Cannot use optimizing backend with less than 2 draws")
  }
  structure(list(model = model, args = args, n_retries_hessian = n_retries_hessian), class = "SBC_backend_rstan_optimizing")
}


#' @export
SBC_fit.SBC_backend_rstan_optimizing <- function(backend, generated, cores) {
  for(attempt in 1:backend$n_retries_hessian) {
    start <- proc.time()
    fit <- do.call(rstan::optimizing,
                   combine_args(list(object = backend$model,
                                     data = generated),
                                backend$args)
    )
    end <- proc.time()
    fit$time <- (end - start)["elapsed"]

    if(fit$return_code != 0) {
      stop("Optimizing was not succesful")
    }
    # This signals production of draws was OK
    if(nrow(fit$theta_tilde) > 1) {
      break;
    }
  }

  fit$n_attempts <- attempt

  if(nrow(fit$theta_tilde) == 1) {
    stop("Optimizing did not return draws.\n",
         "This is most likely due to numerical problems with the Hessian, check model output.\n",
         "You may also consider increasing `n_retries_hessian`")
  }


  structure(fit, class = "RStanOptimizingFit")
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_rstan_optimizing <- function(backend) {
  rlang::hash(list(model = backend$model@model_code, args = backend$args))
}

#' @export
SBC_fit_to_draws_matrix.RStanOptimizingFit <- function(fit) {
  posterior::as_draws_matrix(fit$theta_tilde)
}


#' @export
SBC_backend_iid_draws.SBC_backend_rstan_optimizing <- function(backend) {
  TRUE
}

#' @export
SBC_fit_to_diagnostics.RStanOptimizingFit <- function(fit, fit_output, fit_messages, fit_warnings) {
  res <- data.frame(
    time = fit$time,
    n_attempts = fit$n_attempts
  )

  class(res) <- c("SBC_RStanOptimizing_diagnostics", class(res))
  res
}

#' @export
summary.SBC_RStanOptimizing_diagnostics <- function(x) {
  summ <- list(
    n_fits = nrow(x),
    max_time = max(x$time),
    n_multiple_attempts = sum(x$n_attempts > 1)
  )

  structure(summ, class = "SBC_RStanOptimizing_diagnostics_summary")
}

#' @export
get_diagnostic_messages.SBC_RStanOptimizing_diagnostics <- function(x) {
  get_diagnostic_messages(summary(x))
}


#' @export
get_diagnostic_messages.SBC_RStanOptimizing_diagnostics_summary <- function(x) {
  SBC_diagnostic_messages(
    rbind(
      data.frame(ok = TRUE, message = paste0("Maximum time was ", x$max_time, " sec.")),
      data.frame(ok = x$n_multiple_attempts == 0,
                 message = paste0( x$n_multiple_attempts, " (", round(100 * x$n_multiple_attempts / x$n_fits),
                                   "%) of fits required multiple attempts to produce usable Hessian."))
    )
  )
}
