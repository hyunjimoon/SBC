
#' Plot Rank Histogram given rank array "ranks"
#'
#' @param ranks array of dimension(n_iter, n_pars) where n_iter=number of posterior draw iterations, n_pars the number of parameters of interest
#' @param par names of parameter to plot
#' @param bins number of histogram bins to plot default is 20
#' @import ggplot2
#' @export
plot_hist <- function(ranks, par, bins=20){
  CI = stats::qbinom(c(0.05,0.5,0.95), size=dim(ranks)[1], prob = 1/(bins))
  ggplot() + aes(ranks[, par]) + geom_histogram(bins=bins) +
    geom_hline(yintercept = CI, color="black", linetype="dashed") +
    xlab("rank") + ylab("count") + ggtitle(paste("Rank Histogram for parameter", par))
}

#' Plot ECDF given rank array "ranks"
#' \href{https://arxiv.org/abs/1903.08008}{arxiv::1903.08008} by A. Vehtari et al.
#'
#' @param ranks array of dimension(n_iter, n_pars) where n_iter=number of posterior draw iterations, n_pars the number of parameters of interest
#' @param par names of parameter to plot
#' @import ggplot2
#' @export
plot_ecdf <- function(ranks, par){

  S <- dim(ranks)[1]
  r.scale <- rank(ranks[, par], ties.method="average")
  q95 <- stats::qbeta(0.95, r.scale+1, S - r.scale + 1)
  q05 <- stats::qbeta(0.05, r.scale+1, S - r.scale + 1)
  ggplot() + aes(r.scale / S) + stat_ecdf() +
    geom_line(aes(y=q95), color="blue") + geom_line(aes(y=q05), color="blue") +
    xlab("Fractional Rank") + ylab("ECDF") + ggtitle(paste("ECDF for parameter", par))

}

# plot.ecdf.diff <- function(ranks, par){
#   # Plot ECDF values centered around given rank array "ranks"
#   # Please refer to A. Vehtari et al.
#   #
#   # ranks: array of dimension(n_iter, n_pars) where n_iter=number of posterior draw iterations, n_pars the number of parameters of interest
#   # par: names of parameter to plot
#   S <- dim(ranks)[1]
#   r.scale <- rank(ranks[, par], ties.method="average")  # rescale ranks
#   q95 <- qbeta(0.95, r.scale+1, S - r.scale + 1)
#   q05 <- qbeta(0.05, r.scale+1, S - r.scale + 1)
#
#   ecdf.plot <- ggplot() + aes(x = r.scale / S) + stat_ecdf(n=S, pad=TRUE)  # dummy plot to extract ecdf values
#   ecdf.data <- data.frame(y=layer_data(ecdf.plot)$y, x=layer_data(ecdf.plot)$x)
#   ecdf.data[1, ] <- 0
#   ecdf.data[length(ecdf.data$y), ] <- 1
#
#   plot.data <- data.frame(x=r.scale/S, q95=q95 - (r.scale/S), q05=q05 - (r.scale/S))
#
#   ggplot(plot.data) + aes(x=x) + geom_line(data=ecdf.data, aes(x=x, y=y-x), color="black") + #ylim(-1.0, 1.0) +
#     geom_line(aes(y=q95), color="blue") + geom_line(aes(y=q05), color="blue") + xlab("Fractional Rank") +
#     ylab("ECDF") + ggtitle(paste("Centered ECDF for parameter", par))
#
# }
