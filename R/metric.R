#' Distance between binned draws (rank for SBC) and discrete uniform
#'
#' @param ranks array of dimension (n_iter, n_pars) where n_iter=number of posterior draw iterations, n_pars the number of parameters of interest
#' @param par names of parameter to plot
#' @param thin integer in which thinning was applied
#' @param bins number of bins to use for summary
#' @importFrom stats chisq.test integrate
#' @export
rank2unif <- function(results, par, bins = 20){
  pval <- function(bin_count){
    return(chisq.test(bin_count)$p.value)
  }
  unif_samp <- rep(1/bins, bins)
  binned_samp <- rep(NA, bins)
  for (p in par){
    ranks <- results$stats$rank[names(results$stats$rank) %in% p]
    ranks <- array(ranks, dim=c(length(ranks), 1))
    colnames(ranks) <- c(p)
    S <- results$stats$max_rank[1]
    bin_size <- S / bins
    for (n in 1:dim(ranks)[1]) {
      bin <- ceiling(ranks[n, p] / bin_size)
      binned_samp[bin] <- binned_samp[bin] + 1
    }
    rank_list <- list()
    rank_list[[p]] <- list( max_diff =  round(max_diff(binned_samp, unif_samp),3), wasserstein = round(wasserstein(binned_samp, unif_samp),3)) #pval = round(pval(binned_samp),3),
  }
  return(list(par = par, rank_list = rank_list))
}

#' Summarize relational property of overall prior and posterior draws
#'
#' @param priors A posterior::draws_rvars of dimension(n_iterations=1, n_chains=n_sbc_iterations, n_variables=n_variables) which stores prior draws
#' @param posteriors A posterior::draws_Rvars of dimension(n_iterations=n_posterior_draws, n_chains=n_sbc_iterations, n_variables=n_variables), which stores fitted posterior draws
#' @param par names of parameter to summarize
#' @param bins number of bins for prior and post density
#' @export
# TODO debug wasserstein, max difference
set2set <- function(priors, posteriors, par, bins = 20){
  priors <- draws_of(priors[[par]])
  posteriors <- draws_of(posteriors[[par]])
  breaks <- seq(min(priors, posteriors), max(priors, posteriors), length.out = bins)
  h1 <- hist(priors, breaks = breaks, plot = FALSE)
  h2 <- hist(posteriors,  breaks = breaks, plot = FALSE)
  return(list(CJS = cjs_dist(priors, posteriors), KLD= HistogramTools::kl.divergence(h1, h2), JeffDivg = HistogramTools::jeffrey.divergence(h1, h2),
              Minkowski_2 = HistogramTools::minkowski.dist(h1, h2, 2))) #, max_diff = max_diff(h1$density, h2$density), wass = wasserstein(h1$density, h2$density)
}

##' Max difference between binned samples with the same length
##'
##' @param x numeric vector of density from first distribution
##' @param y numeric vector of density from second distribution
##' @param ... unused
##' @return distance value based on max difference
##' @export
max_diff <- function(x, y){
  diff <- abs(x- y)
  return(diff[which.max(diff)])
}

##' wasserstein distance between binned samples
##'
##' @param x numeric vector of density from first distribution
##' @param y numeric vector of density from second distribution
##' @param ... unused
##' @return distance value based on max difference
##' @export
# TODO need testing
wasserstein <- function(x, y){
  tempf <- Vectorize(function(i) abs((x[i]/sum(x)  - y[i]/sum(y)))) # expected sums = 1
  val <- integrate(tempf,1,5, rel.tol=.Machine$double.eps^.05)$value
  return(val)
}
# wasserstein <- function(x, y, bin_count){
#   bins <- length(bin_count)
#   unif <- rep(1/bins, bins)
#   M <- sum(bin_count)
#   tempf <- Vectorize(function(i)  abs(bin_count[i]/M  - unif[i]))
#   val <- integrate(tempf,1,bins, rel.tol=.Machine$double.eps^.05)$value
#   return(val)
# }

##' Cumulative Jensen-Shannon divergence
##'
##' Computes the cumulative Jensen-Shannon distance between two
##' samples.
##'
##' The Cumulative Jensen-Shannon distance is a symmetric metric based
##' on the cumulative Jensen-Shannon divergence. The divergence CJS(P || Q) between
##' two cumulative distribution functions P and Q is defined as:
##'
##' \deqn{CJS(P || Q) = \sum P(x) \log \frac{P(x)}{0.5 (P(x) + Q(x))} + \frac{1}{2 \ln 2} \sum (Q(x) - P(x))}
##'
##' The symmetric metric is defined as:
##'
##' \deqn{CJS_{dist}(P || Q) = \sqrt{CJS(P || Q) + CJS(Q || P)}}
##'
##' This has an upper bound of \eqn{\sqrt \sum (P(x) + Q(x))}
##'
##' @param x numeric vector of draws from first distribution
##' @param y numeric vector of draws from second distribution
##' @param x_weights numeric vector of weights of first distribution
##' @param y_weights numeric vector of weights of second distribution
##' @param ... unused
##' @return distance value based on CJS computation.
##' @references Nguyen H-V., Vreeken J. (2015).  Non-parametric
##'   Jensen-Shannon Divergence.  In: Appice A., Rodrigues P., Santos
##'   Costa V., Gama J., Jorge A., Soares C. (eds) Machine Learning
##'   and Knowledge Discovery in Databases.  ECML PKDD 2015. Lecture
##'   Notes in Computer Science, vol 9285.  Springer, Cham.
##'   \code{doi:10.1007/978-3-319-23525-7_11}
##' @export
cjs_dist <- function(x, y, x_weights, y_weights, ...) {
  if (class(x)[1] == "rvar"){
    x <- c(draws_of(x))
    #cat("y at cjs")
    #print(y)
    y <- c(draws_of(y))
    x_weights <-  rep(1/length(x), length(x))
    y_weights <-  rep(1/length(y), length(y))
  }
  x_weights <-  rep(1/length(x), length(x))
  y_weights <-  rep(1/length(y), length(y))
  # sort draws and weights
  x_idx <- order(x)
  x <- x[x_idx]
  wp <- x_weights[x_idx]

  y_idx <- order(y)
  y <- y[y_idx]
  wq <- y_weights[y_idx]

  # add end point of final step
  x_v <- x[length(x)] + x[length(x)] - x[length(x) - 1]
  y_v <- y[length(y)] + y[length(y)] - y[length(y) - 1]

  # calculate widths of each step
  x_diff <- diff(c(x, x_v))
  y_diff <- diff(c(y, y_v))

  # calculate required weighted ecdfs
  Px <- spatstat.geom::ewcdf(x, weights = wp)(x)
  Qx <- spatstat.geom::ewcdf(y, weights = wq)(x)

  Py <- spatstat.geom::ewcdf(x, weights = wp)(y)
  Qy <- spatstat.geom::ewcdf(y, weights = wq)(y)

  # calculate integral of ecdfs
  Px_int <- drop(Px %*% x_diff)
  Qx_int <- drop(Qx %*% x_diff)

  Py_int <- drop(Py %*% y_diff)
  Qy_int <- drop(Qy %*% y_diff)

  # calculate cjs
  cjs_PQ <-  x_diff %*% (
    Px * (log(Px, base = 2) -
            log(0.5 * Px + 0.5 * Qx, base = 2)
    )
  ) + 0.5 / log(2) * (Qx_int - Px_int)

  cjs_QP <- y_diff %*% (
    Qy * (log(Qy, base = 2) -
            log(0.5 * Qy + 0.5 * Py, base = 2)
    )
  ) + 0.5 / log(2) * (Py_int - Qy_int)

  # calculate upper bound
  bound <- Px_int + Qy_int

  # normalise with respect to upper bound
  out <- (sqrt(cjs_PQ + cjs_QP)) / sqrt(bound)

  return(drop(out))
}
