
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

# #' Plot ECDF given rank array "ranks"
# #' \href{https://arxiv.org/abs/1903.08008}{arxiv::1903.08008} by A. Vehtari et al.
# #'
# #' @param ranks array of dimension(n_iter, n_pars) where n_iter=number of posterior draw iterations, n_pars the number of parameters of interest
# #' @param par names of parameter to plot
# #' @import ggplot2
# #' @export
# plot_ecdf <- function(ranks, par){
#
#   S <- dim(ranks)[1]
#   r.scale <- rank(ranks[, par], ties.method="average")
#   q95 <- stats::qbeta(0.95, r.scale+1, S - r.scale + 1)
#   q05 <- stats::qbeta(0.05, r.scale+1, S - r.scale + 1)
#   ggplot() + aes(r.scale / S) + stat_ecdf() +
#     geom_line(aes(y=q95), color="blue") + geom_line(aes(y=q05), color="blue") +
#     xlab("Fractional Rank") + ylab("ECDF") + ggtitle(paste("ECDF for parameter", par))
#
# }

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

#' @export
#' @rdname ECDF-plots
#' @param pit For 'ppc_ecdf_intervals' and 'ppc_ecdf_intervals_difference', the
#'   PIT values of one or more samples can be provided directly causing'y' and
#'   'yrep' to be ignored.
#' @param K For 'ppc_ecdf_intervals' and 'ppc_ecdf_intervals_difference',
#'   number of uniformly spaced evaluation points for the ECDF or ECDFs. Affects
#'   the granularity of the plot and can significantly speed up the computation
#'   of the simultaneous confidence bands. Defaults to 'length(y)',
#'   'K = ncol(yrep)', or lastly to 'K = ncol(pit)' depending on which one is
#'   provided.
plot_ecdf <- function(
  sbc_workflow_obj,
  gamma,
  K,
  var,
  ...,
  prob = 0.95,
  size = 1,
  alpha = 0.33
) {

  if(is.null(sbc_workflow_obj$calculated_ranks)){
    stop("No rank data is available. Please run SBCWorkflow$calculate_rank first.")
  }
  if(missing(var)){
    stop("Please specify the parameter name to plot as argument var")
  }

  if(posterior::nvariables(posterior::as_draws_array(posterior::subset_draws(sbc_workflow_obj$prior_samples, var)))){

  }

  pit <- ranks_to_empirical_pit(sbc_workflow_obj$calculated_ranks, posterior::niterations(sbc_workflow_obj$posterior_samples))
  N <- nrow(pit)
  if (missing(K)) {
    K <- N
  }
  if (missing(gamma)) {
    gamma <- adjust_gamma(
      N = N,
      L = 1,
      K = K,
      conf_level = prob
    )
  }
  limits <- ecdf_intervals(
    N = N,
    L = 1,
    K = K,
    gamma = gamma)
  z <- seq(0,1, length.out = K + 1)
  data = data.frame(apply(pit, 2, function(col) ecdf(col)(z) ))
  # construct figure
  ggplot(data) +
    geom_ribbon(
      data = data.frame(limits),
      aes_(
        x = c(0, rep(z[2:(K + 1)], each = 2)),
        ymax = ~ upper / N,
        ymin = ~ lower / N,
        color = "theoretical CDF",
        fill = "theoretical CDF"
      ),
      alpha = alpha,
      size = size) +
    geom_step(
      data = data,
      aes_(x = z, y = data[[var]], color = "sample ECDF"),
      size = size
    ) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_color_manual(
      name = "",
      values = rlang::set_names(
        c("skyblue1", "black"),
        c("theoretical CDF", "sample ECDF")),
      labels = c(
        "theoretical CDF" = expression(italic("theoretical CDF")),
        "sample ECDF" = expression(italic("sample ECDF"))
      )
    ) +
    scale_fill_manual(
      name = "",
      values = c("theoretical CDF" = "skyblue",
                 "sample ECDF" = NA),
      labels = c(
        "theoretical CDF" = expression(italic("theoretical CDF")),
        "sample ECDF" = expression(italic("sample ECDF"))
      )
    ) +
    xlab(NULL) +
    ylab(NULL)
}

#' @export
#' @rdname ECDF-plots
plot_ecdf_diff <- function(
  sbc_workflow_obj,
  gamma,
  K,
  var,
  prob = 0.95,
  size = 1,
  alpha = 0.33
) {
  if(is.null(sbc_workflow_obj$calculated_ranks)){
    stop("No rank data is available. Please run SBCWorkflow$calculate_rank first.")
  }
  pit <- ranks_to_empirical_pit(sbc_workflow_obj$calculated_ranks, posterior::niterations(sbc_workflow_obj$posterior_samples))
  N <- nrow(pit)
  if (missing(K)) {
    K <- N
  }
  if (missing(gamma)) {
    gamma <- adjust_gamma(
      N = N,
      L = 1,
      K = K,
      conf_level = prob
    )
  }
  limits <- ecdf_intervals(
    N = N,
    L = 1,
    K = K,
    gamma = gamma)
  z <- seq(0,1, length.out = K + 1)
  data = data.frame(apply(pit, 2, function(col) ecdf(col)(z) - z))
  ggplot(data) +
    geom_ribbon(
      data = data.frame(limits),
      aes_(
        x = c(0, rep(z[2:(K + 1)], each = 2)),
        ymax = ~ upper / N - c(rep(z[1:K], each = 2), 1),
        ymin = ~ lower / N - c(rep(z[1:K], each = 2), 1),
        color = "theoretical CDF",
        fill = "theoretical CDF"
      ),
      alpha = alpha,
      size = size) +
    geom_step(
      data = data,
      aes_(x = z, y = data[[var]], color = "sample ECDF")
    ) +
    scale_color_manual(
      name = "",
      values = rlang::set_names(
        c("skyblue1", "black"),
        c("theoretical CDF", "sample ECDF")),
      labels = c(
        "theoretical CDF" = expression(italic("theoretical CDF")),
        "sample ECDF" = expression(italic("sample ECDF"))
      )
    ) +
    scale_fill_manual(
      name = "",
      values = c("theoretical CDF" = "skyblue",
                 "sample ECDF" = NA),
      labels = c(
        "theoretical CDF" = expression(italic("theoretical CDF")),
        "sample ECDF" = expression(italic("sample ECDF"))
      )
    ) +
    xlab(NULL) +
    ylab(NULL)
}
