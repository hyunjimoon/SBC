brier_score <- function(x, y) {
  sum((x-y)^2)
}

#' @title Binary calibration tests
#'
#' @param x the predicted success probabilities
#' @param y the actual observed outcomes (just 0 or 1)
#' @param alpha the type I error rate for the test
#' @param B number of boostrap samples for the null distribution
#'
#' @description Dimitriadis et al. propose several tests based on
#' comparing actual predictions to predictions when the probabilities are
#' calibrated. This yields several possible tests of correctly calibrated
#' predictions (i.e. that the expected proportion of true values matches the
#' predicted probability).
#'
#' @details
#' The `brier_` functions represent a test based on brier score, while
#' the `miscalibration_` functions represent a test based on miscalibration.
#' In both cases we evaluate the null distribution via bootstrapping.
#'
#' @returns `brier_resampling_test` and `miscalibration_resampling_test` return
#' an object of class `htest`, `brier_resampling_p` and `miscalibration_resampling_p`
#' return just the p-value (for easier use with automated workflows).
#' `binary_miscalibration` computes just the miscalibration component using
#' the PAV (pool adjacent violators) algorithm.
#'
#'
#' @references     T. Dimitriadis, T. Gneiting, & A.I. Jordan,
#' Stable reliability diagrams for probabilistic classifiers,
#' Proc. Natl. Acad. Sci. U.S.A. 118 (8) e2016191118,
#' https://doi.org/10.1073/pnas.2016191118 (2021).
#'
#' @rdname binary-calibration-tests
#' @export
brier_resampling_test <- function(x, y, alpha = 0.05, B = 10000) {
  dname <- paste0("x = ", deparse1(substitute(x)), ", y = ", deparse1(substitute(y)))

  actual_brier <- brier_score(x, y)
  brier_null <- replicate(B, {
    yrep <- rbinom(length(x), size = 1, prob = x)
    brier_score(x, yrep)
  })

  p <- max(mean(actual_brier <= brier_null), 0.5/B)

  param <- quantile(brier_null, probs = 1 - alpha)
  names(param) <- paste0(scales::percent(1 - alpha), " rejection limit")

  structure(list(
    method = paste0("Bootstrapped binary Brier score test (using ", B, " samples)"),
    data.name = dname,
    p.value = p,
    estimate = c("Brier score" = actual_brier),
    parameter = param
  ),
  class = "htest")
}

#' @rdname binary-calibration-tests
#' @export
brier_resampling_p <- function(x, y, B = 10000) {
  actual_brier <- brier_score(x, y)
  brier_null <- replicate(B, {
    yrep <- rbinom(length(x), size = 1, prob = x)
    brier_score(x, yrep)
  })
  max(mean(actual_brier <= brier_null), 0.5/B)
}

#' @rdname binary-calibration-tests
#' @export
binary_miscalibration <- function(x,y) {
  require_package_version("monotone", "0.1.2", "miscalibration computations")
  ord <- order(x, -y)
  x <- x[ord]
  y <- y[ord]
  #CEP_pav <- stats::isoreg(y)$yf
  CEP_pav <- monotone::monotone(y)
  #Using brier score
  Sc <- mean((CEP_pav - y)^2)
  mean((x - y) ^2) - Sc
}

# Faster reimplementation from https://www.pnas.org/doi/full/10.1073/pnas.2016191118#sec-4
# and the reliabilitydiag package
miscalibration_resampling_nulldist <- function(x,y, B = 1000) {
  replicate(B, {
    yrep <- rbinom(length(x), size = 1, prob = x)
    binary_miscalibration(x, yrep)
  })
}

#' @rdname binary-calibration-tests
#' @export
miscalibration_resampling_p <- function(x,y, B = 10000) {
  actual_miscalibration <- binary_miscalibration(x,y)
  misc_null <- miscalibration_resampling_nulldist(x, y, B)
  max(mean(actual_miscalibration <= misc_null), 0.5/B)
}

#' @rdname binary-calibration-tests
#' @export
miscalibration_resampling_test <- function(x, y, alpha = 0.05, B = 10000) {
  dname <- paste0("x = ", deparse1(substitute(x)), ", y = ", deparse1(substitute(y)))

  actual_miscalibration <- binary_miscalibration(x,y)
  misc_null <- miscalibration_resampling_nulldist(x, y, B)
  p <- max(mean(actual_miscalibration <= misc_null), 0.5/B)

  param <- quantile(misc_null, probs = 1 - alpha)
  names(param) <- paste0(scales::percent(1 - alpha), " rejection limit")

  structure(list(
    method = paste0("Bootstrapped binary miscalibration test (using ", B, " samples)"),
    data.name = dname,
    p.value = p,
    estimate = c("miscalibration" = actual_miscalibration),
    parameter = param
  ),
  class = "htest")
}

gaffke_m <- function(probs, B = 10000) {
  require_package_version("MCMCpack", "1.0.0", "the Gaffke test")
  u_diff <- MCMCpack::rdirichlet(B, alpha = rep(1, length(probs) + 1))

  probs_sort <- sort(probs)
  z_upr <- c(probs_sort, 1)
  m_matrix_upr <- sweep(u_diff, MARGIN = 2, STATS = z_upr, FUN = "*")
  m_upr <- rowSums(m_matrix_upr)

  #stopifnot(identical(sort(1 - probs), rev(1 - probs_sort)))
  z_lwr <- c(rev(1 - probs_sort), 1)
  m_matrix_lwr <- sweep(u_diff, MARGIN = 2, STATS = z_lwr, FUN = "*")
  m_lwr <- rowSums(m_matrix_lwr)

  list(lwr = m_lwr, upr = m_upr)
}

gaffke_ci_from_m <- function(m, alpha = 0.05) {
  m_lwr <- m$lwr
  m_upr <- m$upr

  as.numeric(c(
    1 - quantile(m_lwr, probs = 1 - alpha / 2),
    quantile(m_upr, probs = 1 - alpha / 2)
  ))
}

#' @rdname gaffke_test
#' @export
gaffke_ci <- function(probs, B = 10000, alpha = 0.05) {
  m <- gaffke_m(probs, B, alpha)
  gaffke_ci_from_m(m, alpha)
}

gaffke_p_from_m <- function(m, mu, B, alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)

  m_lwr <- m$lwr
  m_upr <- m$upr

  prob_low <- mean(1-m_lwr <= mu)
  if(prob_low == 0) {
    prob_low <- 0.5/B
  }
  prob_high <- mean(m_upr >= mu)
  if(prob_high == 0) {
    prob_high <- 0.5/B
  }
  if(alternative == "two.sided") {
    return(min(prob_low, prob_high, 0.5) * 2)
  } else if(alternative == "less") {
    return(prob_high)
  } else if(alternative == "greater") {
    return(prob_low)
  } else {
    stop("Invalid alternative")
  }
}

#' @rdname gaffke_test
#' @export
gaffke_p <- function(probs, mu = 0.5, alpha = 0.05, B = 10000, alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)

  m <- gaffke_m(probs, B, alpha)
  gaffke_p_from_m(m, mu, B, alternative)
}

#' Non-parametric test for the mean of a bounded variable.
#'
#' @param x a vector of observed values
#' @param mu the mean under null hypothesis
#' @param alpha the level of the test
#' @param lb the lower bound for `x`
#' @param ub the upper bound for `x`
#' @param B number of bootstrap samples for the null distribution
#' @param alternative the alternative for the test.
#'
#' @details The test is expected to be valid for any bounded distribution without further
#' assumptions. The test has been proven valid only for special cases but
#' no counterexample is known despite some efforts in the literature to find
#' some.
#'
#' @description Test a null hypothesis about the mean of i.i.d. samples.
#' The test is based on Gaffke 2005, though a more detailed analysis and
#' exposition can be found in Learned-Miller & Thomas 2020.
#'
#' @returns `gaffke_test` returns an object of class `htest`, `gaffke_p` and
#' `gaffke_ci` return just the p-value / CI as numeric for easier use in batch
#' workflows.
#'
#' @references Gaffke, N. (2005).
#' “Three test statistics for a nonparametric one-sided hypothesis on the mean
#' of a nonnegative variable.” Mathematical Methods of Statistics, 14(4): 451–467.
#'
#' Learned-Miller, E. and Thomas, P. S. (2020).
#' “A New Confidence Interval for the Mean of a Bounded Random Variable.”
#' https://arxiv.org/abs/1905.06208
#'
#' @rdname gaffke_test
#' @export
gaffke_test <- function(x, mu = 0.5, alpha = 0.05, lb = 0, ub = 1, B = 10000, alternative = c("two.sided", "less", "greater")) {
  dname <- deparse1(substitute(x))
  alternative <- match.arg(alternative)

  stopifnot(length(lb) == 1)
  stopifnot(length(ub) == 1)
  stopifnot(is.finite(lb))
  stopifnot(is.finite(ub))
  stopifnot(all(x >= lb))
  stopifnot(all(x <= ub))
  stopifnot(length(B) == 1 && B > 1)
  stopifnot(0 < alpha && alpha < 1)
  stopifnot(mu >= lb && mu <= ub)

  x_scaled <- (x - lb) / (ub - lb)
  mu_scaled <- (mu - lb) / (ub - lb)
  m <- gaffke_m(x_scaled, B = B)
  p <- gaffke_p_from_m(m, mu_scaled, alternative = alternative)
  ci <- gaffke_ci_from_m(m, alpha = alpha)
  attr(ci, "conf.level") <- 1 - alpha

  structure(list(
    method = paste0("Gaffke's test for the mean of a bounded variable  (using ", B, " samples)"),
    data.name = dname,
    p.value = p,
    alternative = alternative,
    null.value = c("mean" = mu),
    conf.int = ci,
    estimate = c("mean" = mean(x)),
    parameter = c("lower bound" = lb, "upper bound" = ub)
  ),
  class = "htest")
}
