library(cmdstanr)
library(ggplot2)
# prior sample: N samples, for each iteration
# posterior sample: N * Y samples, Y sample in each iteration
# rank = sum(prior < posterior)

# sbc.cmdstanr_posterior2arr(cmdstanfit, pars){
#   draws <- cmdstanfit$draws(variables=pars)
# }


# Basic class to keep data identifiable
#
# @prior: A named array of dimensions(n_iter, n_pars) where n_iter=number of SBC draws, n_pars the number of parameters. Names should be applied as parameter names.
# @posterior: A named array of dimensions(n_iter, n_pars) where n_iter=number of SBC draws, n_pars the number of parameters. Names should be applied as parameter names.
#   Equivalent to posterior.as_draws_array
# @model.name: An user-defined string to help identify what this data came from
SBCData <- setClass("SBCData", 
                    slots = list(prior="array", posterior="array", model.name="character"))


# setMethod("initialize", "SBCData", function(.Object, prior, posterior, model.name){
# })


sbc.rank <- function(SBCData.object, thin){
  # Given prior and posterior samples, generate rank count for each bin.
  #
  # n_iter: number of iterations
  # n_sample: number of posterior samples per iteration
  # n_pars: number of parameters of interest
  # prior array dimensions: (n_iter, n_pars)
  # posterior array dimensions: (n_iter, n_pars, n_sample)
  # return array dimensions: (n_iter, n_pars)
  prior_dims = dim(SBCData.object@prior)
  posterior_dims = dim(SBCData.object@posterior)

  if (prior_dims[1] != posterior_dims[1] || prior_dims[2] != posterior_dims[2]){
    stop("Dimension mismatch error!\n dim(SBCData.object@prior)[1] == dim(SBCData.object@posterior)[1] && dim(SBCData.object@prior)[2] == dim(SBCData.object@posterior)[2] must be satisfied")
  }
  n_iter <- posterior_dims[1]
  n_pars <- posterior_dims[2]
  n_sample <- posterior_dims[3]
  
  par_names <- unlist(dimnames(SBCData.object@prior)[2])
  
  thinner <- seq(from=1, to=n_sample, by=thin)
  ranks <- array(rep(0, n_iter * n_pars), dim=c(n_iter, n_pars))
  dimnames(ranks)[2] <- list(par_names)
  for(i in 1:n_iter){
    for(j in 1:n_pars){
      ranks[i, par_names[j]] <- sum(SBCData.object@prior[i, par_names[j]] < SBCData.object@posterior[i, par_names[j], ][thinner])
    }
  }
  return(ranks)
}

sbc.summary <- function(ranks, par, thin, bins = 20){
  # Summarize Ranks on different distance maeasure metrics
  #
  # ranks: array of dimension(n_iter, n_pars) where n_iter=number of posterior draw iterations, n_pars the number of parameters of interest
  # par: names of parameter to plot
  # thin: integer in which thinning was applied
  # bins=20: number of bins to use for summary
  pval <- function(bin_count){
    return(chisq.test(bin_count)$p.value)
  }
  max_diff <- function(bin_count){
    bins <- length(bin_count)
    diff <- abs(mean(bin_count) - bin_count)
    val <- diff[which.max(diff)] / mean(bin_count)
    return(val)
  }
  wasserstein <- function(bin_count){
    bins <- length(bin_count)
    unif <- rep(1/bins, bins)
    M <- sum(bin_count)
    tempf <- Vectorize(function(i)  abs(bin_count[i]/M  - unif[i]))
    val <- integrate(tempf,1,bins, rel.tol=.Machine$double.eps^.05)$value
    return(val)
  }
  S <- dim(ranks)[1]
  bin_size <- S / bins
  bin_count <- rep(0, bins)
  for (s in 1:S) {
    bin <- ceiling(ranks[s, par] / bins)
    bin_count[bin] <- bin_count[bin] + 1
  }
  print(paste0("pval: ", round(pval(bin_count),3), " max_diff: ", round(max_diff(bin_count),3), " wasserstein: ", round(wasserstein(bin_count),3)))
}

sbc.plot.hist <- function(ranks, par, thin, bins=20){
  # Plot Rank Histogram given rank array "ranks"
  #
  # ranks: array of dimension(n_iter, n_pars) where n_iter=number of posterior draw iterations, n_pars the number of parameters of interest
  # par: names of parameter to plot
  # thin: integer in which thinning was applied
  # bins=20: number of histogram bins to plot
  CI = qbinom(c(0.05,0.5,0.95), size=dim(ranks)[1], prob = 1/(bins))
  ggplot() + aes(ranks[, par]) + geom_histogram(bins=bins) + geom_hline(yintercept = CI, color="black", linetype="dashed") +
    xlab("rank") + ylab("count") + ggtitle(paste("Rank Histogram for parameter", par))
}

sbc.plot.ecdf <- function(ranks, par){
  # Plot ECDF given rank array "ranks"
  # https://arxiv.org/abs/1903.08008 by A. Vehtari et al.
  #
  # ranks: array of dimension(n_iter, n_pars) where n_iter=number of posterior draw iterations, n_pars the number of parameters of interest
  # par: names of parameter to plot
  S <- dim(ranks)[1]
  r.scale <- rank(ranks[, par], ties.method="average")
  q95 <- qbeta(0.95, r.scale+1, S - r.scale + 1)
  q05 <- qbeta(0.05, r.scale+1, S - r.scale + 1)
  ggplot() + aes(r.scale / S) + stat_ecdf() +
    geom_line(aes(y=q95), color="blue") + geom_line(aes(y=q05), color="blue") + xlab("Fractional Rank") + ylab("ECDF") + ggtitle(paste("ECDF for parameter", par))
  
}

sbc.plot.ecdf.diff <- function(ranks, par){
  # Plot ECDF values centered around given rank array "ranks"
  # Please refer to A. Vehtari et al.
  #
  # ranks: array of dimension(n_iter, n_pars) where n_iter=number of posterior draw iterations, n_pars the number of parameters of interest
  # par: names of parameter to plot
  S <- dim(ranks)[1]
  r.scale <- rank(ranks[, par], ties.method="average")  # rescale ranks
  q95 <- qbeta(0.95, r.scale+1, S - r.scale + 1)
  q05 <- qbeta(0.05, r.scale+1, S - r.scale + 1)

  ecdf.plot <- ggplot() + aes(x = r.scale / S) + stat_ecdf(n=S, pad=TRUE)  # dummy plot to extract ecdf values
  ecdf.data <- data.frame(y=layer_data(ecdf.plot)$y, x=layer_data(ecdf.plot)$x)
  ecdf.data[1, ] <- 0
  ecdf.data[length(ecdf.data$y), ] <- 1

  plot.data <- data.frame(x=r.scale/S, q95=q95 - (r.scale/S), q05=q05 - (r.scale/S))
  
  ggplot(plot.data) + aes(x=x) + geom_line(data=ecdf.data, aes(x=x, y=y-x), color="black") + #ylim(-1.0, 1.0) + 
    geom_line(aes(y=q95), color="blue") + geom_line(aes(y=q05), color="blue") + xlab("Fractional Rank") +
    ylab("ECDF") + ggtitle(paste("Centered ECDF for parameter", par))
  
}

