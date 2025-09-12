SBC_nuts_diagnostic_types <- function() {
  list(
    max_chain_time = SBC_numeric_diagnostic("maximum time per chain", report = "max"),
    n_divergent = SBC_count_diagnostic("divergences", report = "max", lower_thresh = 0),
    n_max_treedepth = SBC_count_diagnostic("iterations that saturated max treedepth", lower_thresh = 0, label_short = "max treedepths"),
    n_rejects = SBC_count_diagnostic("steps rejected", report = "max", lower_thresh = 0),
    min_bfmi = SBC_numeric_diagnostic("E-BFMI", report = "min", upper_thresh = 0.2),
    n_failed_chains = SBC_count_diagnostic("failed chains", lower_thresh = 0)
  )
}


#' @export
summary.SBC_nuts_diagnostics <- function(diagnostics) {
  summ <- list(
    n_fits = nrow(diagnostics),
    max_chain_time = max(diagnostics$max_chain_time),
    has_divergent = sum(diagnostics$n_divergent > 0),
    max_divergent = max(diagnostics$n_divergent),
    has_treedepth = sum(diagnostics$n_max_treedepth > 0),
    max_max_treedepth = max(diagnostics$n_max_treedepth),
    has_rejects = sum(diagnostics$n_rejects > 0),
    max_rejects = max(diagnostics$n_rejects)
  )

  if(!is.null(diagnostics$min_bfmi)) {
    summ$has_low_bfmi = sum(is.na(diagnostics$min_bfmi) | diagnostics$min_bfmi < 0.2)
  }

  if(!is.null(diagnostics$n_failed_chains)) {
    if(any(is.na(diagnostics$n_failed_chains))) {
      problematic_sim_ids <- paste0(which(is.na(diagnostics$n_failed_chains)), collapse = ", ")
      warning("Fits for simulations ", problematic_sim_ids, " had NA for n_failed_chains.")
    }
    summ$has_failed_chains = sum(is.na(diagnostics$n_failed_chains) | diagnostics$n_failed_chains > 0)
  }

  structure(summ, class = "SBC_nuts_diagnostics_summary")
}


get_diagnostic_messages.SBC_nuts_diagnostics <- function(x) {
  get_diagnostic_messages(summary(x))
}

# Not currently used. Discussion at https://discourse.mc-stan.org/t/summarising-rhat-values-over-multiple-variables-fits/23957/
# was helpful and need to figure out a better way to handle this.
#
# Use the Generalized extreme value distribution
# to get a quantile of maximum of `n_vars` random values
# distributed as N(1, 0.005).
# Following https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution#Example_for_Normally_distributed_variables
# The approximation looks good for n_vars >= 10
# for smaller, we just plug in a log function with appropriate scale
# See https://math.stackexchange.com/questions/89030/expectation-of-the-maximum-of-gaussian-random-variables
# for a discussion on why log
get_expected_max_rhat <- function(n_vars, prob = 0.99, approx_sd = 0.005) {
  stopifnot(is.numeric(n_vars))
  stopifnot(all(n_vars >= 1))

  # Maximum of n_vars standardized normals
  gumbel_approx <- function(n) {
    # Gumbel params
    mu_n <- qnorm(1 - 1/n)
    sigma_n <- qnorm(1 - (1 / n) * exp(-1)) - mu_n

    # Inverse CDF of gumbel with xi = 0
    mu_n - (sigma_n * log(-log(prob)))
  }


  linear_bound <- 10
  approx_at_bound <- gumbel_approx(linear_bound)
  value_at_1 <- qnorm(prob)
  linear_scale <- (approx_at_bound - value_at_1 ) / log(linear_bound)

  std_val_max <- dplyr::if_else(n_vars < linear_bound,
                                value_at_1 + linear_scale * log(n_vars),
                                gumbel_approx(n_vars)
  )
  1 + std_val_max * approx_sd
}

#' @export
get_diagnostic_messages.SBC_nuts_diagnostics_summary <- function(x) {
  message_list <- list()
  i <- 1
  if(!is.null(x$has_failed_chains)) {
    if(x$has_failed_chains > 0) {
      msg <- paste0(x$has_failed_chains, " (", round(100 * x$has_failed_chains / x$n_fits), "%) fits had some failed chains.")
      message_list[[i]] <- data.frame(ok = FALSE, message = msg)
    } else {
      message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had failed chains.")
    }
    i <- i + 1
  }

  if(x$has_divergent > 0) {
    msg <- paste0(x$has_divergent, " (", round(100 * x$has_divergent / x$n_fits),
                  "%) fits had divergent transitions. Maximum number of divergences was ",
                  x$max_divergent, ".")
    message_list[[i]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had divergent transitions.")
  }
  i <- i + 1

  if(x$has_treedepth > 0) {
    msg <- paste0(x$has_treedepth, " (", round(100 * x$has_treedepth / x$n_fits),
                  "%) fits had iterations that saturated max treedepth. Maximum number of max treedepth was ",
                  x$max_max_treedepth, ".")
    message_list[[i]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had iterations that saturated max treedepth.")
  }
  i <- i + 1

  if(!is.null(x$has_low_bfmi)) {
    if(x$has_low_bfmi > 0) {
      msg <- paste0(x$has_low_bfmi, " (", round(100 * x$has_low_bfmi / x$n_fits), "%) fits had low BFMI.")
      message_list[[i]] <- data.frame(ok = FALSE, message = msg)
    } else {
      message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had low BFMI.")
    }
    i <- i + 1
  }

  if(x$has_rejects > 0) {
    msg <- paste0(x$has_rejects, " (", round(100 * x$has_rejects / x$n_fits), "%) fits had some steps ",
                  "rejected. Maximum number of rejections was ", x$max_rejects, ".")
    message_list[[i]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[i]] <- data.frame(ok = TRUE, message = "No fits had steps rejected.")
  }
  i <- i + 1

  message_list[[i]] <- data.frame(ok = TRUE, message = paste0("Maximum time per chain was ", x$max_chain_time, " sec."))
  i <- i + 1

  SBC_diagnostic_messages(do.call(rbind, message_list))
}

#' @export
print.SBC_nuts_diagnostics_summary <- function(x) {
  msg <- get_diagnostic_messages(x)
  print(msg)
  invisible(x)
}


#' @export
SBC_ADVI_diagnostics_types <- function() {
  list(
    elbo_converged = SBC_logical_diagnostic("ELBO converged", error_value = FALSE),
    n_rejects = SBC_count_diagnostic("some steps rejected", report = "max", lower_thresh = 0),
    time = SBC_numeric_diagnostic("time", report = "max")
  )
}


#' @export
summary.SBC_ADVI_diagnostics <- function(x) {
  summ <- list(
    n_fits = nrow(x),
    n_elbo_not_converged = sum(!x$elbo_converged),
    has_rejects = sum(x$n_rejects > 0),
    max_rejects = max(x$n_rejects),
    max_time = max(x$time)

  )

  structure(summ, class = "SBC_ADVI_diagnostics_summary")
}

#' @export
get_diagnostic_messages.SBC_ADVI_diagnostics <- function(x) {
  get_diagnostic_messages(summary(x))
}


#' @export
get_diagnostic_messages.SBC_ADVI_diagnostics_summary <- function(x) {
  message_list <- list()

  if(x$n_elbo_not_converged == 0) {
    message_list[[1]] <-
      data.frame(ok = TRUE, message = "All fits converged.")
  } else {
    message_list[[1]] <-
      data.frame(ok = FALSE,
                 message = paste0(
                   x$n_elbo_not_converged, " (", round(100 * x$n_elbo_not_converged / x$n_fits),
                   "%) of fits did not converge."))
  }


  if(x$has_rejects > 0) {
    msg <- paste0(x$has_rejects, " (", round(100 * x$has_rejects / x$n_fits), "%) fits had some steps ",
                  "rejected. Maximum number of rejections was ", x$max_rejects, ".")
    message_list[[2]] <- data.frame(ok = FALSE, message = msg)
  } else {
    message_list[[2]] <- data.frame(ok = TRUE, message = "No fits had steps rejected.")
  }

  message_list[[3]] <- data.frame(ok = TRUE, message = paste0("Maximum time was ", x$max_time, " sec."))

  SBC_diagnostic_messages(do.call(rbind, message_list))
}
