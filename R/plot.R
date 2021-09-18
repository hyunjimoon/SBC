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
#' @import ggplot2
plot_rank_hist.data.frame <- function(x, parameters = NULL, bins = NULL, prob = 0.95, max_rank = x$max_rank) {
  if(!all(c("parameter", "rank") %in% names(x))) {
    stop("The data.frame needs a 'parameter' and 'rank' columns")
  }
  n_simulations <- dplyr::summarise(dplyr::group_by(x, parameter), count = dplyr::n())$count
  if(length(unique(n_simulations)) > 1) {
    stop("Differing number of SBC steps per parameter not supported.")
  }

  if(is.null(max_rank)) {
    stop("max_rank either has to be supplied explicitly or be a column in the data")
  }
  max_rank <- unique(max_rank)
  if(length(max_rank) > 1) {
    stop("Differing max_rank across parameters is not supported yet.")
  }

  n_simulations <- unique(n_simulations)

  if(is.null(bins)){
    bins <- guess_bins(max_rank, n_simulations)
  } else if(bins > max_rank + 1) {
    stop("Cannot use more bins than max_rank + 1")
  }

  if(!is.null(parameters)) {
    x <- dplyr::filter(x, parameter %in% parameters)
  }

  if(nrow(x) == 0) {
    stop("No data for the selected parameters.")
  }

  #CI - taken from https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R


  # Bins can differ by size (at most by 1). Build a CI that is conservative,
  # i.e. includes lower quantile of smalelr bins and higher quantile of larger bins
  larger_bin_size <- ceiling((max_rank / bins))
  smaller_bin_size <- floor((max_rank / bins))
  CI = qbinom(c(0.5 * (1 - prob),0.5,0.5 * (1 + prob)), size=n_simulations,prob  =  larger_bin_size / max_rank)
  ci_lower = qbinom(0.5 * (1 - prob), size=n_simulations,prob  =  smaller_bin_size / max_rank)
  ci_mean = qbinom(0.5, size=n_simulations,prob  =  1 / bins)
  ci_upper = qbinom(0.5 * (1 + prob), size=n_simulations,prob  =  larger_bin_size / max_rank)

  CI_polygon_x <- c(-0.1*max_rank,0,-0.1*max_rank,1.1 * max_rank,max_rank,1.1 * max_rank,-0.1 * max_rank)
  CI_polygon_y <- c(ci_lower,ci_mean,ci_upper,ci_upper,ci_mean,ci_lower,ci_lower)

  #The visualisation style taken as well from   https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R
  ggplot(x, aes(x = rank)) +
          geom_segment(aes(x=0,y=ci_mean,xend=max_rank,yend=ci_mean),colour="grey25") +
          geom_polygon(data=data.frame(x= CI_polygon_x,y= CI_polygon_y),aes(x=x,y=y),fill="skyblue",color="skyblue1",alpha=0.33) +
          geom_histogram(breaks =  seq(0, max_rank, length.out = bins + 1), closed = "left" ,fill="#808080",colour="black") +
          scale_y_continuous("count") +
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
  min(max_rank + 1, max(floor(N / 10), 5))
}

#' Plot the ECDF-based plots.
#'
#'
#' See the methods for [data_for_ecdf_plots()] for available data formats.
#'
#' \href{https://arxiv.org/abs/1903.08008}{arxiv::1903.08008} by A. Vehtari et al.
#' @export
#' @rdname ECDF-plots
#' @param x object supporting the [data_for_ecdf_plots()] method.
#' @param parameters optional subset of parameters to show in the plot
#' @param gamma TODO
#' @param prob the width of the plotted confidence interval for the ECDF.
#' @param size size passed to [ggplot2::geom_ribbon()] for the confidence band
#' @param alpha alpha level of the confidence band
#' @param K number of uniformly spaced evaluation points for the ECDF or ECDFs. Affects
#'   the granularity of the plot and can significantly speed up the computation
#'   of the simultaneous confidence bands. Defaults to the smaller of number of
#'   ranks per parameter and the maximum rank.
#' @param ... additional arguments passed to [data_for_ecdf_plots()].
#' Most notably, if `x` is matrix, a `max_rank` parameter needs to be given.
#' @import ggplot2
#' @seealso [plot_coverage()]
plot_ecdf <- function(x,
                      parameters = NULL,
                      K = NULL,
                      gamma = NULL,
                      prob = 0.95,
                      size = 1,
                      alpha = 0.33, ...) {

  ecdf_data <-
    data_for_ecdf_plots(x, parameters = parameters,
                        prob = prob, K = K, gamma = gamma, ...)

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
#' @import ggplot2
plot_ecdf_diff <- function(x,
                           parameters = NULL,
                           K = NULL,
                           gamma = NULL,
                           prob = 0.95,
                           size = 1,
                           alpha = 0.33, ...) {
  ecdf_data <-
    data_for_ecdf_plots(x, parameters = parameters,
                        prob = prob, K = K, gamma = gamma, ...)

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

  if(is.null(max_rank)) {
    stop("max_rank either has to be supplied explicitly or be a column in the data")
  }
  max_rank <- unique(max_rank)
  if(length(max_rank) > 1) {
    stop("Differing max_rank across parameters is not supported yet.")
  }

  summ <- dplyr::summarise(dplyr::group_by(stats, parameter), count = dplyr::n(), .groups = "drop")
  if(length(unique(summ$count)) > 1) {
    stop("Not all variables have the same number of simulations.")
  }

  rank <- dplyr::select(stats, dataset_id, parameter, rank)
  rank_matrix <- tidyr::pivot_wider(rank, names_from = "parameter",
                                    values_from = "rank")
  rank_matrix <- as.matrix(dplyr::select(rank_matrix, -dataset_id))


  data_for_ecdf_plots(rank_matrix, max_rank = max_rank, prob = prob,
                      gamma = gamma, K = K)
}

data_for_ecdf_plots.matrix <- function(x,
                                               max_rank,
                                       parameters = NULL,
                                                 prob = 0.95,
                                                 gamma = NULL,
                                                 K = NULL,
                                                 size = 1,
                                                 alpha = 0.33) {
  ranks_matrix <- x
  if(any(!is.finite(ranks_matrix))) {
    stop("Ranks may only contain finite values")
  }

  if(!is.null(parameters)) {
    ranks_matrix <- ranks_matrix[, parameters]
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


#' Prior/posterior contraction plot.
#'
#' The rationale for this plot and its interpretaion is explained in
#' Mike Betancourt's
#' [Towards A Principled Bayesian Workflow](https://betanalpha.github.io/assets/case_studies/principled_bayesian_workflow.html#132_A_Bayesian_Eye_Chart).
#'
#' @param x object containing results (a data.frame or [SBC_results()] object).
#' @param prior_sd a named vector of prior standard deviations for your parameters.
#' Either pass in analytically obtained values or use [calculate_prior_sd()] to get an empirical estimate from
#' an `SBC_datasets` object.
#' @param parameters parameters to show in the plot or `NULL` to show all
#' must correspond a field already computed in the results (most likely `"mean"` and `"median"`).
#' @param scale which scale of variability you want to see - either `"sd"` for standard deviation
#' or `"var"` for variance.
#' @param alpha the alpha for the points
#' @return a ggplot2 plot object
#' @export
plot_contraction <- function(x, prior_sd, parameters = NULL, scale = "sd", alpha = 0.8) {
  UseMethod("plot_contraction")
}

#' @export
plot_contraction.SBC_results <- function(x, prior_sd, parameters = NULL, scale = "sd", alpha = 0.8) {
  plot_contraction(x$stats, prior_sd = prior_sd, parameters = parameters, alpha = alpha)
}

#' @export
plot_contraction.data.frame <- function(x, prior_sd, parameters = NULL, scale = "sd", alpha = 0.8) {
  if(!all(c("parameter", "sd") %in% names(x))) {
    stop("The data.frame needs a 'parameter' and 'sd' columns")
  }

  if(!is.numeric(prior_sd) || is.null(names(prior_sd))) {
    stop("prior_sd has to be a named vector")
  }

  if(!is.null(parameters)) {
    prior_sd <- prior_sd[names(prior_sd) %in% parameters]
    x <- dplyr::filter(x, parameter %in% parameters)
  }

  if(nrow(x) == 0 || length(prior_sd) == 0) {
    stop("No data to plot.")
  }

  shared_params <- intersect(unique(x$parameter), names(prior_sd))
  if(length(shared_params) < length(unique(x$parameter))) {
    warning("Some parameters do not have prior_sd in the data: ", setdiff(unique(x$parameter), shared_params))
  }
  if(length(shared_params) < length(prior_sd)) {
    warning("Some prior_sd values do not have counterpart in the data: ", setdiff(names(prior_sd), shared_params))
  }

  x <- dplyr::filter(x, parameter %in% shared_params)

  x$prior_sd <- prior_sd[x$parameter]
  if(scale == "sd") {
    x <- dplyr::mutate(x, contraction = 1 - sd / prior_sd)
  } else if(scale == "var") {
    x <- dplyr::mutate(x, contraction = 1 - (sd / prior_sd)^2)
  }

  ggplot2::ggplot(x, aes(x = contraction, y = z_score)) + geom_point(alpha = alpha) +
    expand_limits(x = c(0,1)) +
    facet_wrap(~parameter)
}


#' Plot the simulated "true" values versus posterior estimates
#'
#' @param x object containing results (a data.frame or [SBC_results()] object).
#' @param parameters parameters to show in the plot or `NULL` to show all
#' @param estimate which estimate to use for the central tendency,
#' must correspond a field already computed in the results (most likely `"mean"` and `"median"`).
#' @param uncertainty which estimates to use for uncertainty (a character vector of length 2)
#' must correspond a field already computed in the results. Pass `NULL` to avoid showing uncertainty at all.
#' @param alpha the alpha for the points and uncertainty intervals
#' @return a ggplot2 plot object
#' @export
plot_sim_estimated <- function(x, parameters = NULL, estimate = "mean",
                               uncertainty = c("q5", "q95"),
                               alpha = 0.8) {
  UseMethod("plot_sim_estimated")
}

#' @export
plot_sim_estimated.SBC_results <- function(x, parameters = NULL, estimate = "mean",
                                           uncertainty = c("q5", "q95"),
                                           alpha = 0.8) {
  plot_sim_estimated(x$stats, parameters = parameters, estimate = estimate,
                     uncertainty = uncertainty, alpha = alpha)
}

#' @export
plot_sim_estimated.data.frame <- function(x, parameters = NULL, estimate = "mean",
                                          uncertainty = c("q5", "q95"),
                                          alpha = 0.8) {
  if(!all(c("parameter", estimate, uncertainty) %in% names(x))) {
    stop("The data.frame needs a 'parameter' and '", estimate, "' columns")
  }

  if(!is.null(parameters)) {
    x <- dplyr::filter(x, parameter %in% parameters)
  }

  x$estimate__ <- x[[estimate]]

  if(!is.null(uncertainty)) {
    if(length(uncertainty) != 2) {
      stop("'uncertainty' has to be null or a character vector of length 2")
    }
    x$low__ <- x[[uncertainty[1]]]
    x$high__ <- x[[uncertainty[2]]]
    all_aes <- aes(x = simulated_value, y = estimate__, ymin = low__, ymax = high__)
    main_geom <- geom_pointrange(alpha = alpha, fatten = 1.5)
  } else {
    main_geom <- geom_point(alpha = alpha)
    all_aes <- aes(x = simulated_value, y = estimate__)
  }

  if(nrow(x) == 0) {
    stop("No data to plot.")
  }

  ggplot2::ggplot(x, all_aes) +
    geom_abline(intercept = 0, slope = 1, color = "skyblue1", size = 2) +
    main_geom +
    scale_y_continuous(estimate) +
    facet_wrap(~parameter, scales = "free")
}


#' Plot the observed coverage and its uncertainty
#'
#' Please refer to [empirical_coverage()] for details on computation
#' and limitations of this plot as well as details on the arguments.
#'
#' @param x object containing results (a data.frame or [SBC_results()] object).
#' @param parameters parameters to show in the plot or `NULL` to show all
#' @param prob the with of the uncertainty interval to be shown
#' @return a ggplot2 plot object
#' @export
plot_coverage <- function(x, parameters = NULL, prob = 0.95,
                          interval_type = "central") {
  UseMethod("plot_coverage")
}

#' @rdname plot_coverage
#' @export
plot_coverage.SBC_results <- function(x, parameters = NULL, prob = 0.95,
                                      interval_type = "central") {
  plot_coverage(x$stats, parameters = parameters, prob = prob, interval_type = interval_type)
}

#' @rdname plot_coverage
#' @export
plot_coverage.data.frame <- function(x, parameters = NULL, prob = 0.95,
                                     interval_type = "central") {
  if(!all(c("parameter", "rank", "max_rank") %in% names(x))) {
    stop(SBC_error("SBC_invalid_argument_error",
                   "The stats data.frame needs a 'parameter', 'rank' and 'max_rank' columns"))
  }

  if(!is.null(parameters)) {
    x <- dplyr::filter(x, parameter %in% parameters)
  }

  max_max_rank <- max(x$max_rank)
  coverage <- empirical_coverage(x, (0:max_max_rank) / (max_max_rank + 1), prob = prob,
                                interval_type = interval_type)

  ggplot2::ggplot(coverage, aes(x = width_represented, y = estimate,
                                ymin = ci_low, ymax = ci_high)) +
    geom_ribbon(fill = "black", alpha = 0.33) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, color = "skyblue1", size = 2) +
    #geom_abline(intercept = 0, slope = 1, color = "skyblue1", size = 2) +
    geom_line() +
    scale_x_continuous(paste0(interval_type, " interval width")) +
    scale_y_continuous("Observed coverage") +
    facet_wrap(~parameter)


}
