#' @export
diagnostic_types.SBC_nuts_diagnostics <- function(diags) {
  types <- list(
    max_chain_time = numeric_diagnostic("time per chain", report = "max", digits = 1),
    n_divergent = count_diagnostic("divergences", error_above = 0),
    n_max_treedepth = count_diagnostic("iterations that saturated max treedepth", error_above = 0, label_short = "max treedepths"),
    n_rejects = count_diagnostic("steps rejected", error_above = 0),
    min_bfmi = numeric_diagnostic("E-BFMI", report = "min", error_below = 0.2, digits = 3)
  )

  if("n_failed_chains" %in% names(diags)) {
    types$n_failed_chains <- count_diagnostic("failed chains", error_above = 0)
  }

  types
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
diagnostic_types.SBC_ADVI_diagnostics <- function(diags) {
  list(
    elbo_converged = logical_diagnostic("ELBO converged", error_value = FALSE),
    n_rejects = count_diagnostic("some steps rejected", error_above = 0),
    time = numeric_diagnostic("time", report = "max")
  )
}

