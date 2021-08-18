
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

#' Plot rank histogram of an SBC results.
#'
#' By default the support is for `SBC_results` objects and data frames in the same
#' format as the `$stats` element of `SBC_results`.
#' @param x Object supporting the plotting method.
#' @export
plot_rank_hist <- function(x, parameters = NULL, bins = NULL, prob = 0.95, ...) {
  UseMethod("plot_rank_hist")
}

#' @export
plot_rank_hist.data.frame <- function(x, parameters = NULL, bins = NULL, prob = 0.95, max_rank = x$max_rank) {
  if(!all(c("parameter", "rank") %in% names(x))) {
    stop("The data.frame needs a 'parameter' and 'rank' columns")
  }
  n_simulations <- dplyr::summarise(dplyr::group_by(x, parameter), count = n())$count
  if(length(unique(n_simulations)) > 1) {
    stop("Differing number of SBC steps per parameter not supported.")
  }

  max_rank <- unique(max_rank)
  if(length(max_rank) > 1) {
    stop("Differing max_rank across parameters is not supported yet.")
  }

  n_simulations <- unique(n_simulations)

  if(is.null(bins)){
    bins <- guess_bins(max_rank, n_simulations)
  }

  if(!is.null(parameters)) {
    x <- dplyr::filter(x, parameter %in% parameters)
  }

  if(nrow(x) == 0) {
    stop("No data for the selected parameters.")
  }

  #CI - taken from https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R


  larger_bin_size <- ceiling((max_rank / bins))
  CI = qbinom(c(0.5 * (1 - prob),0.5,0.5 * (1 + prob)), size=n_simulations,prob  =  larger_bin_size / max_rank)
  ci_lower = CI[1]
  ci_mean = CI[2]
  ci_upper = CI[3]

  CI_polygon_x <- c(-0.1*max_rank,0,-0.1*max_rank,1.1 * max_rank,max_rank,1.1 * max_rank,-0.1 * max_rank)
  CI_polygon_y <- c(ci_lower,ci_mean,ci_upper,ci_upper,ci_mean,ci_lower,ci_lower)

  #The visualisation style taken as well from   https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R
  x %>%
          ggplot(aes(x = rank)) +
          geom_segment(aes(x=0,y=ci_mean,xend=max_rank,yend=ci_mean),colour="grey25") +
          geom_polygon(data=data.frame(x= CI_polygon_x,y= CI_polygon_y),aes(x=x,y=y),fill="skyblue",color="skyblue1",alpha=0.33) +
          geom_histogram(breaks =  seq(0, max_rank, length.out = bins + 1), closed = "left" ,fill="#808080",colour="black") +
          facet_wrap(~parameter, scales = "free_y")

}


#' @export
plot_rank_hist.SBC_results <- function(x, parameters = NULL, bins = NULL, prob = 0.95) {
  x <- validate_SBC_results(x)
  max_rank <- unique(x$stats$max_rank)
  if(length(max_rank) > 1) {
    stop("Differing max_rank across parameters not supported yet.")
  }

  plot_rank_hist(x$stats, parameters = parameters, bins = bins, prob = prob, max_rank = max_rank)
}

#' Guess the number of bins for [plot_rank_hist()].
#' @param N the number of ranks observed
#' @param max_rank the maximum rank observed
guess_bins <- function(max_rank, N) {
  min(max_rank, max(floor(N / 10), 5))
}

#' Plot the ECDF-based plots.
#'
#'
#' \href{https://arxiv.org/abs/1903.08008}{arxiv::1903.08008} by A. Vehtari et al.
#' @export
#' @rdname ECDF-plots
#' @param x object supporting the [data_for_ecdf_plots()] method.
#' @param parameters optional subset of parameters to show in the plot
#' @param gamma TODO
#' @param prob the width of the plotted confidence interval for the ECDF.
#' @param size size passed to [ggplot::geom_ribbon()] for the confidence band
#' @param alpha alpha level of the confidence band
#' @param K number of uniformly spaced evaluation points for the ECDF or ECDFs. Affects
#'   the granularity of the plot and can significantly speed up the computation
#'   of the simultaneous confidence bands. Defaults to the smaller of number of
#'   ranks per parameter and the maximum rank.
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
  data_for_ecdf_plots(x$stats, parameters = parameters, prob = prob, gamma = gamma, K = K)
}


data_for_ecdf_plots.data.frame <- function(x, parameters = NULL,
                                           prob = 0.95,
                                           gamma = NULL,
                                           K = NULL,
                                           max_rank = x$max_rank) {
  stats <- x
  if(!is.null(parameters)) {
    stats <- dplyr::filter(stats, parameter %in% parameters)
  }

  max_rank <- unique(max_rank)
  if(length(max_rank) > 1) {
    stop("Differing max_rank across parameters is not supported yet.")
  }

  summ <- dplyr::summarise(dplyr::group_by(stats, parameter), count = dplyr::n(), .groups = "drop")
  if(length(unique(summ$count)) > 1) {
    stop("Not all varaibles have the same number of ranks")
  }

  rank <- dplyr::select(stats, dataset_id, parameter, rank)
  rank_matrix <- tidyr::pivot_wider(rank, names_from = "parameter",
                                    values_from = "rank")
  rank_matrix <- as.matrix(dplyr::select(rank_matrix, -dataset_id))


  data_for_ecdf_plots(rank_matrix, max_rank = max_rank, prob = prob,
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
