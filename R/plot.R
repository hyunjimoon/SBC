
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
plot_ecdf <- function(x,
                      parameters = NULL,
                      K = NULL,
                      gamma = NULL,
                      prob = 0.95,
                      size = 1,
                      alpha = 0.33) {

  ecdf_data <-
    data_for_ecdf_plots(x, parameters = parameters,
                        prob = prob, K = K, gamma = gamma)

  N <- ecdf_data$N
  K <- ecdf_data$K
  z <- ecdf_data$z

  ecdf_df <- dplyr::mutate(ecdf_data$ecdf_df, type = "sample ECDF")
  limits_df <- ecdf_data$limits_df
  limits_df_trans <- data.frame(
    x = c(0, rep(z[2:(K + 1)], each = 2)),
    ymax = limits_df$upper / N,
    ymin = limits_df$lower / N,
    type = "theoretical CDF"
  )

  # construct figure
  ggplot(ecdf_df) +
    geom_ribbon(
      data = limits_df_trans,
      aes(x = x, ymax = ymax, ymin = ymin, color = type, fill = type),
      alpha = alpha,
      size = size) +
    geom_step(
      aes(x = z, y = ecdf, color = type)
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
    ylab(NULL) +
    facet_wrap(~ parameter)
}

#' @export
#' @rdname ECDF-plots
plot_ecdf_diff <- function(x,
                           parameters = NULL,
                           K = NULL,
                           gamma = NULL,
                           prob = 0.95,
                           size = 1,
                           alpha = 0.33) {
  ecdf_data <-
    data_for_ecdf_plots(x, parameters = parameters,
                        prob = prob, K = K, gamma = gamma)

  N <- ecdf_data$N
  K <- ecdf_data$K
  z <- ecdf_data$z

  ecdf_df <- dplyr::mutate(ecdf_data$ecdf_df, z_diff = ecdf - z, type = "sample ECDF")
  limits_df <- ecdf_data$limits_df
  limits_df_trans <- data.frame(
    x = c(0, rep(z[2:(K + 1)], each = 2)),
    ymax = limits_df$upper / N - c(rep(z[1:K], each = 2), 1),
    ymin = limits_df$lower / N - c(rep(z[1:K], each = 2), 1),
    type = "theoretical CDF"
  )
  ggplot(ecdf_df) +
    geom_ribbon(
      data = limits_df_trans,
      aes(x = x, ymax = ymax, ymin = ymin, color = type, fill = type),
      alpha = alpha,
      size = size) +
    geom_step(
      aes(x = z, y = z_diff, color = type)
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
    ylab(NULL) +
    facet_wrap(~ parameter)
}



#' Maybe not export in the end? Useful for debugging
#' @export
data_for_ecdf_plots <- function(x, ...,
                                        prob = 0.95,
                                        gamma = NULL,
                                        K = NULL
                                ) {
  UseMethod("data_for_ecdf_plots")
}


data_for_ecdf_plots.SBC_results <- function(x, parameters = NULL,
                                                    prob = 0.95,
                                                    gamma = NULL,
                                                    K = NULL) {
  stats <- x$stats
  if(!is.null(parameters)) {
    stats <- dplyr::filter(stats, parameter %in% parameters)
  }

  if(length(unique(stats$max_rank)) > 1) {
    stop("Not all parameters have the same max_rank")
  }

  summ <- dplyr::summarise(dplyr::group_by(stats, parameter), count = dplyr::n(), .groups = "drop")
  if(length(unique(summ$count)) > 1) {
    stop("Not all varaibles have the same number of ranks")
  }

  rank <- dplyr::select(stats, dataset_id, parameter, rank)
  rank_matrix <- tidyr::pivot_wider(rank, names_from = "parameter",
                                              values_from = "rank")
  rank_matrix <- as.matrix(dplyr::select(rank_matrix, -dataset_id))


  data_for_ecdf_plots(rank_matrix, max_rank = stats$max_rank[1], prob = prob,
                              gamma = gamma, K = K)
}


data_for_ecdf_plots.SBCWorkflow <- function(
  x, parameters = NULL,
  prob = 0.95,
  gamma = NULL,
  K = NULL
) {
  sbc_workflow_obj <- x
  if(is.null(sbc_workflow_obj$calculated_ranks)){
    stop("No rank data is available. Please run SBCWorkflow$calculate_rank first.")
  }

  data_for_ecdf_plots(sbc_workflow_obj$calculated_ranks,
                      max_rank = posterior::niterations(sbc_workflow_obj$posterior_samples),
                          gamma = gamma, K = K, prob = prob)


}

data_for_ecdf_plots.matrix <- function(x,
                                               max_rank,
                                                 prob = 0.95,
                                                 gamma = NULL,
                                                 K = NULL,
                                                 size = 1,
                                                 alpha = 0.33) {
  ranks_matrix <- x
  if(any(!is.finite(ranks_matrix))) {
    stop("Ranks may only contain finite values")
  }
  pit <- ranks_to_empirical_pit(ranks_matrix, max_rank)
  N <- nrow(pit)
  if (is.null(K)) {
    K <- min(max_rank + 1, N)
  }
  if (is.null(gamma)) {
    gamma <- adjust_gamma(
      N = N,
      L = 1,
      K = K,
      conf_level = prob
    )
  }
  limits_df <- as.data.frame(ecdf_intervals(
    N = N,
    L = 1,
    K = K,
    gamma = gamma))
  z <- seq(0,1, length.out = K + 1)

  ecdf_vals <- apply(pit, 2, function(col) ecdf(col)(z))

  ecdf_df <- as.data.frame(ecdf_vals)
  ecdf_df$..z <- z
  ecdf_df <- tidyr::pivot_longer(ecdf_df, -..z, names_to = "parameter", values_to = "ecdf")
  ecdf_df <- dplyr::rename(ecdf_df, z = ..z)

  structure(list(limits_df = limits_df, ecdf_df = ecdf_df, K = K, N = N, z = z),
            class = "SBC_ecdf_data")
}
