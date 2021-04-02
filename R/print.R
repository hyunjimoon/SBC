
#' Summarize Ranks on different distance maeasure metrics
#'
#' @param ranks array of dimension(n_iter, n_pars) where n_iter=number of posterior draw iterations, n_pars the number of parameters of interest
#' @param par names of parameter to plot
#' @param thin integer in which thinning was applied
#' @param bins number of bins to use for summary
#' @importFrom stats chisq.test integrate
#' @export
print_summary <- function(ranks, par, thin, bins = 20){
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
  print(paste0("pval: ", round(pval(bin_count),10), " max_diff: ", round(max_diff(bin_count),3), " wasserstein: ", round(wasserstein(bin_count),3)))
}
