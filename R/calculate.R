

#' Given prior and posterior samples, generate rank count for each bin.
#'
#' @param prior A named array of dimensions(n_iter, n_pars) where n_iter=number of SBC draws, n_pars the number of parameters. Names should be applied as parameter names.
#' @param posterior A named array of dimensions(n_samples, n_pars, n_iter) where n_iter=number of SBC draws, n_pars the number of parameters, and n_samples the number of samples per SBC draw. Names should be applied as parameter names. Equivalent to posterior.as_draws_array
#' @param thin Integer in which thinning samples will be applied
#' prior array dimensions: (n_iter, n_pars)
#' posterior array dimensions: (n_sample, n_pars, n_iter)
#'
#' @return array of dimensions(n_iter, n_pars)
#' @export
calculate_rank <- function(prior, posterior, thin){

  prior_dims = dim(prior)
  posterior_dims = dim(posterior)

  if (prior_dims[1] != posterior_dims[3]){
    stop(paste("Dimension mismatch error!",
               "dim(prior)[1] == dim(posterior)[3] must be satisfied",
               paste("prior dimensions:", prior_dims[1], prior_dims[2]),
               paste("posterior dimensions:", posterior_dims[1], posterior_dims[2], posterior_dims[3]),
               "", sep="\n"))
  }
  n_sample <- posterior_dims[1]
  n_iter <- posterior_dims[3]

  par_names <- intersect(unlist(dimnames(prior)[2]), unlist(dimnames(posterior)[2]))
  n_pars <- length(par_names)
  if(n_pars == 0){
    stop("There isn't a parameter name both present in the prior and posterior column names. SBC cannot continue. length(intersect(prior, posterior)) == 0")
  }

  thinner <- seq(from=1, to=n_sample, by=thin)
  ranks <- array(rep(0, n_iter * n_pars), dim=c(n_iter, n_pars))
  dimnames(ranks)[2] <- list(par_names)
  for(i in 1:n_iter){
    for(j in 1:n_pars){
      ranks[i, par_names[j]] <- sum(posterior[, par_names[j], i][thinner] < prior[i, par_names[j]])
    }
  }
  return(ranks)
}

adjust_gamma <- function(N, L, K=N, conf_level=0.95) {
  if (any(c(K, N, L) < 1)) {
    abort("Parameters 'N', 'L' and 'K' must be positive integers.")
  }
  if (conf_level >= 1 || conf_level <= 0) {
    abort("Value of 'conf_level' must be in (0,1).")
  }
  if (L==1) {
    gamma <- adjust_gamma_optimize(N, K, conf_level)
  }
  else {
    gamma <- adjust_gamma_simulate(N, L, K, conf_level)
  }
  gamma
}

# Adjust coverage parameter to find silmultaneous confidence intervals for the
# ECDF of a sample from the uniform distribution.
# N - length of samples
# K - number of equally spaced evaluation points, i.e. the right ends of the
# partition intervals.
adjust_gamma_optimize <- function(N, K, conf_level=0.95) {
  target <- function(gamma, conf_level, N, K) {
    z <- 1:(K - 1) / K
    z1 <- c(0,z)
    z2 <- c(z,1)

    # pre-compute quantiles and use symmetry for increased efficiency.
    x2_lower <- qbinom(gamma / 2, N, z2)
    x2_upper <- c(N - rev(x2_lower)[2:K], 1)

    # Compute the total probability of trajectories inside the confidence
    # intervals. Initialize the set and corresponding probasbilities known
    # to be 0 and 1 for the starting value z1 = 0.
    x1 <- 0
    p_int <- 1
    for (i in seq_along(z1)) {
      tmp <- p_interior(
        p_int, x1 = x1, x2 = x2_lower[i]: x2_upper[i],
        z1 = z1[i], z2 = z2[i], gamma = gamma, N = N
      )
      x1 <- tmp$x1
      p_int <- tmp$p_int
    }
    abs(conf_level - sum(p_int))
  }
  optimize(target, c(0, 1 - conf_level), conf_level, N = N, K = K)$minimum
}

# Adjust coverage parameter to find silmultaneous confidence intervals for the
# ECDFs of multiple samples (chains) from the uniform distribution.
# N - length of samples (chains).
# L - number of samples (chains).
# K - number of equally spaced evaluation points, i.e. the right ends of the
# partition intervals.
# M - number of simulations used to determine the 'conf_level' middle quantile.
adjust_gamma_simulate <-function(N, L, K, conf_level=0.95, M=5000) {
  gamma <- numeric(M)
  z <- (1:(K - 1)) / K
  if (L > 1){
    n <- N * (L - 1)
    k <- floor(z * N * L)
    for (m in seq_len(M)) {
      u = u_scale(replicate(L, runif(N)))
      scaled_ecdfs <- apply(outer(u, z, "<="), c(2,3), sum)
      gamma[m] <- 2 * min(
        apply(
          scaled_ecdfs, 1, phyper, m = N, n = n, k = k
        ),
        apply(
          scaled_ecdfs - 1, 1, phyper, m = N, n = n, k = k, lower.tail = FALSE
        )
      )
    }

  }
  else {
    for (m in seq_len(M)) {
      u <- runif(N)
      scaled_ecdf <- colSums(outer(u, z, "<="))
      gamma[m] <- 2 * min(
        pbinom(scaled_ecdf, N, z),
        pbinom(scaled_ecdfs - 1, N, z, lower.tail = FALSE)
      )
    }
  }
  alpha_quantile(gamma, 1 - conf_level)
}

p_interior <- function(p_int, x1, x2, z1, z2, gamma, N) {
  z_tilde <- (z2 - z1) / (1 - z1)

  N_tilde <- rep(N - x1, each = length(x2))
  p_int <- rep(p_int, each = length(x2))
  x_diff <- outer(x2, x1, "-")
  p_x2_int <- p_int * dbinom(x_diff, N_tilde, z_tilde)

  list(p_int = rowSums(p_x2_int), x1 = x2)
}

# 100 * `alpha` percent of the trials are allowed to be rejected.
# In case of ties, return the largest value dominating at most
# 100 * (alpha + tol) percent of the values.
alpha_quantile <- function(gamma, alpha, tol = 0.001) {
  a <- unname(quantile(gamma, probs = alpha))
  a_tol <- unname(quantile(gamma, probs = alpha + tol))
  if (a == a_tol) {
    if (min(gamma) < a) {
      # take the largest value that doesn't exceed the tolerance.
      a <- max(gamma[gamma < a])
    }
  }
  a
}

# Compute simultaneous confidence intervals for one or more samples from the
# standard uniform distribution.
# N - sample length
# L - number of samples
# K - size of uniform partition defining the ECDF evaluation points.
# gamma - coverage parameter for the marginal distribution (binomial for
# one sample and hypergeometric for multiple rank transformed chains).
ecdf_intervals <- function(N, L, K, gamma) {
  lims <- list()
  z <- seq(0,1, length.out = K + 1)
  if (L == 1) {
    lims$lower <- qbinom(gamma / 2, N, z)
    lims$upper <- qbinom(1 - gamma / 2, N, z)
  } else {
    n <- N * (L - 1)
    k <- floor(z * L * N)
    lims$lower <- qhyper(gamma / 2, N, n, k)
    lims$upper <- qhyper(1 - gamma / 2, N, n, k)
  }
  lims$lower <- c(rep(lims$lower[1:K], each=2), lims$lower[K + 1])
  lims$upper <- c(rep(lims$upper[1:K], each=2), lims$upper[K + 1])
  lims
}

# Transform observations in 'x' into their corresponding fractional ranks.
u_scale <- function(x) {
  array(rank(x) / length(x), dim = dim(x), dimnames = dimnames(x))
}

# for each value in 'y', compute the fractional ranks (empirical pit values)
# with respect to 'yrep'.
empirical_pit <- function(y, yrep) {
  (1 +  apply(outer(yrep, y, "<="), 3, sum)) / (1 +length(yrep))
}
