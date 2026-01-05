filter_stats_by_variables_and_hidden <- function(x, variables, show_hidden = FALSE) {
  if(!is.null(variables)) {
    x <- dplyr::filter(x, variable %in% variables)
  } else if("attributes" %in% names(x) & !show_hidden) {
    x <- dplyr::filter(x, !attribute_present_stats(hidden_var_attribute(), attributes))
  }

  if(nrow(x) == 0) {
    stop("No variables to plot.")
  }
  x
}

#' Plot rank histogram of an SBC results.
#'
#' The expected uniform distribution and an approximate confidence interval
#' is also shown. The confidence interval cannot be taken too seriously
#' as it is derived assuming the bins are independent (which they are not).
#' The [plot_ecdf()] and [plot_ecdf_diff()] plots provide better confidence interval
#' but are somewhat less interpretable. See `vignette("rank_visualizations")` for
#' more details.
#'
#' By default the support is for `SBC_results` objects and data frames in the same
#' format as the `$stats` element of `SBC_results`.
#'
#' @param x Object supporting the plotting method.
#' @param variables Names of variables to show
#' @param bins number of bins to be used in the histogram, if left unspecified,
#'   it is determined by [guess_rank_hist_bins()].
#' @param prob The width of the approximate confidence interval shown.
#' @export
plot_rank_hist <- function(x, variables = NULL, bins = NULL, prob = 0.95, ..., parameters = NULL) {
  UseMethod("plot_rank_hist")
}

#' @export
#' @rdname plot_rank_hist
#' @param show_hidden Show variables marked with [hidden_var_attribute()] (by default, those are not shown)
#' @param facet_args extra arguments for the call to [ggplot2::facet_wrap()] when rendering the plot.
#' @import ggplot2
plot_rank_hist.data.frame <- function(x, variables = NULL, bins = NULL, prob = 0.95, max_rank = x$max_rank, show_hidden = FALSE, parameters = NULL, facet_args = list()) {
  # Ensuring backwards compatibility
  if("parameter" %in% names(x)) {
    if(!("variable" %in% names(x))) {
      warning("The x parameter contains a `parameter` column, which is deprecated, use `variable` instead.")
      x$variable <- x$parameter
    }
  }


  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }



  if(!all(c("variable", "rank") %in% names(x))) {
    stop("The data.frame needs a 'variable' and 'rank' columns")
  }

  x <- filter_stats_by_variables_and_hidden(x, variables, show_hidden)

  if(nrow(x) == 0) {
    stop("No data for the selected variables.")
  }

  n_sims <- dplyr::summarise(dplyr::group_by(x, variable), count = dplyr::n())$count
  if(length(unique(n_sims)) > 1) {
    stop("Differing number of SBC steps per variable not supported.")
  }

  if(is.null(max_rank)) {
    stop("max_rank either has to be supplied explicitly or be a column in the data")
  }
  max_rank <- unique(max_rank)
  if(length(max_rank) > 1) {
    stop("Differing max_rank across variables is not supported yet.")
  }

  n_sims <- unique(n_sims)

  if(is.null(bins)){
    bins <- guess_rank_hist_bins(max_rank, n_sims)
  } else if(bins > max_rank + 1) {
    stop("Cannot use more bins than max_rank + 1")
  }

  #CI - taken from https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R


  # Bins can differ by size (at most by 1). Build a CI that is conservative,
  # i.e. includes lower quantile of smalelr bins and higher quantile of larger bins
  larger_bin_size <- ceiling(((max_rank + 1) / bins))
  smaller_bin_size <- floor(((max_rank + 1) / bins))
  ci_lower = qbinom(0.5 * (1 - prob), size=n_sims,prob  =  smaller_bin_size / max_rank)
  ci_mean = qbinom(0.5, size=n_sims,prob  =  1 / bins)
  ci_upper = qbinom(0.5 * (1 + prob), size=n_sims,prob  =  larger_bin_size / max_rank)

  CI_polygon_x <- c(-0.1*max_rank,0,-0.1*max_rank,1.1 * max_rank,max_rank,1.1 * max_rank,-0.1 * max_rank)
  CI_polygon_y <- c(ci_lower,ci_mean,ci_upper,ci_upper,ci_mean,ci_lower,ci_lower)

  all_facet_args <- c(list(~variable, scales = "free_y"), facet_args)

  #The visualisation style taken as well from   https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R
  ggplot(x, aes(x = rank)) +
          geom_segment(aes(x=0,y=ci_mean,xend=max_rank,yend=ci_mean),colour="grey25") +
          geom_polygon(data=data.frame(x= CI_polygon_x,y= CI_polygon_y),aes(x=x,y=y),fill="skyblue",color="skyblue1",alpha=0.33) +
          geom_histogram(breaks =  seq(0, max_rank, length.out = bins + 1), closed = "left" ,fill="#808080",colour="black") +
          labs(y = "count") +
          do.call(facet_wrap, all_facet_args)

}


#' @export
plot_rank_hist.SBC_results <- function(x, variables = NULL, bins = NULL, prob = 0.95, parameters = NULL, ...) {
  x <- validate_SBC_results(x)

  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  max_rank <- unique(x$stats$max_rank)
  if(length(max_rank) > 1) {
    stop("Differing max_rank across variables not supported yet.")
  }

  plot_rank_hist(x$stats, variables = variables, bins = bins, prob = prob, max_rank = max_rank, ...)
}

#' Guess the number of bins for [plot_rank_hist()].
#' @param N the number of ranks observed
#' @param max_rank the maximum rank observed
guess_rank_hist_bins <- function(max_rank, N) {
  min(max_rank + 1, max(floor(N / 10), 5))
}

#' Plot the ECDF-based plots.
#'
#'
#' See `vignette("rank_visualizations")` for
#' more details.
#' See the methods for [data_for_ecdf_plots()] for available data formats.
#'
#' \href{https://arxiv.org/abs/1903.08008}{arxiv::1903.08008} by A. Vehtari et al.
#' @export
#' @rdname ECDF-plots
#' @param x object supporting the [data_for_ecdf_plots()] method.
#' @param variables optional subset of variables to show in the plot
#' @param gamma TODO
#' @param prob the width of the plotted confidence interval for the ECDF.
#' @param size size (linewidth) passed to [ggplot2::geom_ribbon()] for the confidence band
#' @param alpha alpha level of the confidence band
#' @param K number of uniformly spaced evaluation points for the ECDF or ECDFs. Affects
#'   the granularity of the plot and can significantly speed up the computation
#'   of the simultaneous confidence bands. Default value is chosen heuristically.
#'   You can also use `"max"` to represent
#'   the number of ranks or `"min"` to choose a lower but still sensible value.
#' @param combine_variables optionally specify a named list where each entry is a character
#'   vectors which specifies a group of variables that will be displayed in the
#'   same panel. Panel title will be the name of the list element.
#'   A function that takes a character vector as an input and produces a list
#'   can also be specified (see [combine-functions]).
#' @param ecdf_alpha the alpha level of the empirical CDF. Can be either a single number or
#'   a function taking the number of variables that were combined (when `combine_variables`
#'   is specified) and returns a number. By default, plots showing many
#'   ECDFs will have reduced alpha.
#' @param show_hidden Show variables marked with [hidden_var_attribute()]
#'    (by default, those are not shown, available only when `x` is a data.frame)
#' @param facet_args extra arguments for the call to [ggplot2::facet_wrap()] when rendering the plot.
#' @param ... additional arguments passed to [data_for_ecdf_plots()].
#' Most notably, if `x` is matrix, a `max_rank` parameter needs to be given.
#' @param parameters DEPRECATED, use `variables` instead.
#' @import ggplot2
#' @seealso [plot_coverage()]
plot_ecdf <- function(x,
                      variables = NULL,
                      K = NULL,
                      gamma = NULL,
                      prob = 0.95,
                      size = 1,
                      alpha = 0.33,
                      combine_variables = NULL,
                      ecdf_alpha = NULL,
                      show_hidden = FALSE,
                      facet_args = list(),
                      ...,
                      parameters = NULL) {

  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  ecdf_data <-
    data_for_ecdf_plots(x, variables = variables,
                        prob = prob, K = K, gamma = gamma,
                        combine_variables = combine_variables,
                        ecdf_alpha = ecdf_alpha,
                        show_hidden = show_hidden, ...)

  N <- ecdf_data$N
  K <- ecdf_data$K
  z <- ecdf_data$z

  ecdf_df <- dplyr::mutate(ecdf_data$ecdf_df, type = "sample ECDF")
  limits_df <- ecdf_data$limits_df
  limits_df$type <- "theoretical CDF"


  all_facet_args <- c(list(~group), facet_args)

  # construct figure
  ggplot(ecdf_df, aes(color = type, fill = type)) +
    geom_ribbon(
      data = limits_df,
      aes(x = x, ymax = upper, ymin = lower),
      alpha = alpha,
      linewidth = size) +
    geom_step(
      aes(x = z, y = ecdf, group = variable, alpha = alpha)
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
                 "sample ECDF" = "transparent"),
      labels = c(
        "theoretical CDF" = expression(italic("theoretical CDF")),
        "sample ECDF" = expression(italic("sample ECDF"))
      )
    ) +
    scale_alpha_identity() +
    xlab(NULL) +
    ylab(NULL) +
    do.call(facet_wrap, all_facet_args)
}

#' @export
#' @rdname ECDF-plots
#' @import ggplot2
plot_ecdf_diff <- function(x,
                           variables = NULL,
                           K = NULL,
                           gamma = NULL,
                           prob = 0.95,
                           size = 1,
                           alpha = 0.33,
                           combine_variables = NULL,
                           ecdf_alpha = NULL,
                           show_hidden = FALSE,
                           facet_args = list(),
                           ...,
                           parameters = NULL) {
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  ecdf_data <-
    data_for_ecdf_plots(x, variables = variables,
                        prob = prob, K = K, gamma = gamma,
                        combine_variables = combine_variables,
                        ecdf_alpha = ecdf_alpha,
                        show_hidden = show_hidden,
                        ...)

  if(ecdf_data$N < 50 && is.null(K)) {
    message("With less than 50 simulations, we recommend using plot_ecdf as it has better fidelity.\n",
            "Disable this message by explicitly setting the K parameter. ",
            "You can use the strings \"max\" (high fidelity) or \"min\" (nicer plot) or choose a specific integer.")
  }
  N <- ecdf_data$N
  K <- ecdf_data$K
  z <- ecdf_data$z

  ecdf_df <- dplyr::mutate(ecdf_data$ecdf_df, z_diff = ecdf - z, type = "sample ECDF")
  limits_df_trans <- dplyr::mutate(ecdf_data$limits_df,
    ymax = upper - uniform_val,
    ymin = lower - uniform_val,
    type = "theoretical CDF"
  )

  all_facet_args <- c(list(~group, scales = "free_y"), facet_args)

  ggplot(ecdf_df, aes(color = type, fill = type)) +
    geom_ribbon(
      data = limits_df_trans,
      aes(x = x, ymax = ymax, ymin = ymin),
      alpha = alpha,
      linewidth = size) +
    geom_step(
      aes(x = z, y = z_diff, group = variable, alpha = alpha)
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
                 "sample ECDF" = "transparent"),
      labels = c(
        "theoretical CDF" = expression(italic("theoretical CDF")),
        "sample ECDF" = expression(italic("sample ECDF"))
      )
    ) +
    scale_alpha_identity() +
    xlab(NULL) +
    ylab(NULL) +
    do.call(facet_wrap, all_facet_args)
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


#' @export
data_for_ecdf_plots.SBC_results <- function(x, variables = NULL,
                                                    prob = 0.95,
                                                    gamma = NULL,
                                                    K = NULL,
                                                    ...,
                                            parameters = NULL) {
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  data_for_ecdf_plots(x$stats, variables = variables, prob = prob,
  gamma = gamma, K = K, ...)
}


#' @export
data_for_ecdf_plots.data.frame <- function(x, variables = NULL,
                                           prob = 0.95,
                                           gamma = NULL,
                                           K = NULL,
                                           max_rank = x$max_rank,
                                           show_hidden = FALSE,
                                           ...,
                                           parameters = NULL) {

  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  if("parameter" %in% names(x)) {
    if(!("variable" %in% names(x))) {
      warning("The x parameter contains a `parameter` column, which is deprecated, use `variable` instead.")
      x$variable <- x$parameter
    }
  }

  if("dataset_id" %in% names(x)) {
    if(!("sim_id" %in% names(x))) {
      warning("The x parameter contains a `dataset_id` column, which is deprecated, use `sim_id` instead.")
      x$sim_id <- x$dataset_id
    }
  }


  if(!all(c("variable", "rank", "sim_id") %in% names(x))) {
    stop(SBC_error("SBC_invalid_argument_error",
                   "The stats data.frame needs a 'variable', 'rank' and 'sim_id' columns"))
  }

  stats <- filter_stats_by_variables_and_hidden(x, variables, show_hidden)

  if(is.null(max_rank)) {
    stop("max_rank either has to be supplied explicitly or be a column in the data")
  }
  max_rank <- unique(max_rank)
  if(length(max_rank) > 1) {
    stop("Differing max_rank across variables is not supported yet.")
  }

  summ <- dplyr::summarise(dplyr::group_by(stats, variable), count = dplyr::n(), .groups = "drop")
  if(length(unique(summ$count)) > 1) {
    stop("Not all variables have the same number of simulations.")
  }

  rank <- dplyr::select(stats, sim_id, variable, rank)
  rank_matrix <- tidyr::pivot_wider(rank, names_from = "variable",
                                    values_from = "rank")
  rank_matrix <- as.matrix(dplyr::select(rank_matrix, -sim_id))


  data_for_ecdf_plots(rank_matrix, max_rank = max_rank, prob = prob,
                      gamma = gamma, K = K, ...)
}

#' @export
data_for_ecdf_plots.matrix <- function(x,
                                       max_rank,
                                       variables = NULL,
                                       prob = 0.95,
                                       gamma = NULL,
                                       K = NULL,
                                       size = 1,
                                       alpha = 0.33,
                                       combine_variables = NULL,
                                       ecdf_alpha = NULL,
                                       ...,
                                       parameters = NULL) {

  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  ranks_matrix <- x
  if(any(!is.finite(ranks_matrix))) {
    stop("Ranks may only contain finite values")
  }

  if(!is.null(variables)) {
    ranks_matrix <- ranks_matrix[, variables]
  }

  N <- nrow(ranks_matrix)
  if (is.null(K)) {
    if(N < 50) {
      K <- max_rank + 1
    } else {
      K <- min(max_rank + 1, N)
    }
  } else if (K == "max") {
    K <- max_rank + 1
  } else if (K == "min") {
    K <- min(max_rank + 1, N, 100)
  } else if (!is.numeric(K) & length(K) == 1) {
    stop("K must be either a single number, \"max\", \"min\" or NULL")
  } else {
    K <- as.integer(K)
  }
  if (is.null(gamma)) {
    gamma <- adjust_gamma(
      N = N,
      L = 1,
      K = K,
      conf_level = prob
    )
  }
  z <- seq(0,1, length.out = K + 1)
  z_twice <- c(0, rep(z[2:(K + 1)], each = 2))

  limits_df <- as.data.frame(ecdf_intervals(
    N = N,
    L = 1,
    K = K,
    gamma = gamma))
  limits_df <- dplyr::mutate(limits_df,
                             x = z_twice,
                             lower = lower / N,
                             upper = upper / N,
                             # The uniform_val needs to be shifted w.r.t z_twice
                             uniform_val =  c(rep(z[1:K], each = 2), 1))

  # Combining pit and ecdf calculations in one function to avoid
  # numerical problems causing issue #79
  base_vals <- floor((0:K) * ((max_rank + 1) / K))
  ecdf_vals <- matrix(nrow = K + 1, ncol = ncol(ranks_matrix))
  colnames(ecdf_vals) <- colnames(ranks_matrix)
  for(i in 1:(K + 1)) {
    # Note: for pit calculations we would use (col + 1) / (max_rank + 1)
    # For ecdf we would use pit <= base_val / (max_rank + 1)
    # So the "+ 1" and "<=" can be subsumed in "<"
    ecdf_vals[i,] <- colMeans(ranks_matrix < base_vals[i])
  }


  ecdf_df <- as.data.frame(ecdf_vals)
  ecdf_df$..z <- z
  ecdf_df <- tidyr::pivot_longer(ecdf_df, -..z, names_to = "variable", values_to = "ecdf")
  ecdf_df <- dplyr::rename(ecdf_df, z = ..z)
  # Allow user-specified grouping of variables + alpha on ecdf line (issue #88)
  if(is.null(ecdf_alpha)) {
    ecdf_alpha <- \(x) sqrt(1/x)
  } else if(is.numeric(ecdf_alpha) & length(ecdf_alpha) == 1) {
    ecdf_alpha_numeric <- ecdf_alpha
    ecdf_alpha <- \(x) ecdf_alpha_numeric
  } else if(!is.function(ecdf_alpha) | nargs(ecdf_alpha) != 1) {
    stop("`ecdf_alpha` must be a function taking a single argument or a single numerical value")
  }
  if (!is.null(combine_variables)) {
    if(is.function(combine_variables)) {
      combine_variables <- combine_variables(unique(ecdf_df$variable))
    }
    if(!is.list(combine_variables) | is.null(names(combine_variables))) {
      stop("`combine_variables` must be a named list or a function returning a named list")
    }

    if(!identical(unique(table(unlist(combine_variables))), 1L)) {
      stop("Duplicated variable names are not allowed in `combine_variables`")
    }
    if(!all(unlist(combine_variables) %in% ecdf_df$variable)) {
      stop("The following variables in `combine_variables` couldn't be found: ",
        paste(unlist(combine_variables)[!unlist(combine_variables) %in% ecdf_df$variable], collapse = ", "))
    }
    display_names <- names(combine_variables)
    for (i in seq_along(combine_variables)) {
      ecdf_df[ecdf_df$variable %in% combine_variables[[i]], "group"] <- display_names[i]
    }
    ecdf_df$group <- factor(ecdf_df$group, levels = display_names, ordered = TRUE)
    ecdf_df <- dplyr::mutate(ecdf_df,
      alpha = ecdf_alpha(length(unique(variable))), .by = group)
  } else {
    ecdf_df$alpha <- ecdf_alpha(1)
    ecdf_df$group <- ecdf_df$variable
  }

  structure(list(limits_df = limits_df, ecdf_df = ecdf_df, K = K, N = N, z = z_twice),
            class = "SBC_ecdf_data")
}

#' Helper functions to be passed to [ECDF-plots] to combine variables in a single
#' panel.
#'
#' `combine_all_variables` will merge all variables in a single plot, while
#' `combine_array_elements` will merge all elements of any array into a single
#' panel of the plot
#' @param x parameter names
#' @export
#' @rdname combine-functions
combine_all_variables <- function(x) {
  list(all = x)
}

#' @export
#' @rdname combine-functions
combine_array_elements <- function(x) {
  indices_removed <- gsub(r"(\[[^]]*\])", "[]", x)
  unique_arrays <- sort(unique(indices_removed))
  res <- list()
  for(i in 1:length(unique_arrays)) {
    res[[unique_arrays[i]]] <- x[indices_removed == unique_arrays[i]]
  }
  res
}

#' Prior/posterior contraction plot.
#'
#' The rationale for this plot and its interpretation is explained in
#' Mike Betancourt's
#' [Towards A Principled Bayesian Workflow](https://betanalpha.github.io/assets/case_studies/principled_bayesian_workflow.html#132_A_Bayesian_Eye_Chart).
#'
#' @param x object containing results (a data.frame or [SBC_results()] object).
#' @param prior_sd a named vector of prior standard deviations for your variables.
#' Either pass in analytically obtained values or use [calculate_prior_sd()] to get an empirical estimate from
#' an `SBC_datasets` object.
#' @param variables variables to show in the plot or `NULL` to show all
#' must correspond a field already computed in the results (most likely `"mean"` and `"median"`).
#' @param scale which scale of variability you want to see - either `"sd"` for standard deviation
#' or `"var"` for variance.
#' @param alpha the alpha for the points
#' @param parameters DEPRECATED, use `variables` instead.
#' @return a ggplot2 plot object
#' @export
plot_contraction <- function(x, prior_sd, variables = NULL, scale = "sd", alpha = 0.8, ...,  parameters = NULL) {
  UseMethod("plot_contraction")
}

#' @export
plot_contraction.SBC_results <- function(x, prior_sd, variables = NULL, scale = "sd", alpha = 0.8, parameters = NULL) {
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  plot_contraction(x$stats, prior_sd = prior_sd, variables = variables, alpha = alpha)
}

#' @param show_hidden Show variables marked with [hidden_var_attribute()]
#'    (by default, those are not shown, available only when `x` is a data.frame)
#' @export
#' @rdname plot_contraction
plot_contraction.data.frame <- function(x, prior_sd, variables = NULL, scale = "sd", alpha = 0.8, show_hidden = FALSE, parameters = NULL) {
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  # Ensuring backwards compatibility
  if("parameter" %in% names(x)) {
    if(!("variable" %in% names(x))) {
      warning("The x parameter contains a `parameter` column, which is deprecated, use `variable` instead.")
      x$variable <- x$parameter
    }
  }

  if(!all(c("variable", "sd") %in% names(x))) {
    stop("The data.frame needs a 'variable' and 'sd' columns")
  }

  if(!is.numeric(prior_sd) || is.null(names(prior_sd))) {
    stop("prior_sd has to be a named vector")
  }

  x <- filter_stats_by_variables_and_hidden(x, variables, show_hidden)

  if(!is.null(variables)) {
    prior_sd <- prior_sd[names(prior_sd) %in% variables]
  }

  if(nrow(x) == 0 || length(prior_sd) == 0) {
    stop("No data to plot.")
  }

  shared_vars <- intersect(unique(x$variable), names(prior_sd))
  if(length(shared_vars) < length(unique(x$variable))) {
    warning("Some variables do not have prior_sd in the data: ", setdiff(unique(x$variable), shared_vars))
  }
  if(length(shared_vars) < length(prior_sd)) {
    warning("Some prior_sd values do not have counterpart in the data: ", setdiff(names(prior_sd), shared_vars))
  }

  x <- dplyr::filter(x, variable %in% shared_vars)

  x$prior_sd <- prior_sd[x$variable]
  if(scale == "sd") {
    x <- dplyr::mutate(x, contraction = 1 - sd / prior_sd)
  } else if(scale == "var") {
    x <- dplyr::mutate(x, contraction = 1 - (sd / prior_sd)^2)
  }

  ggplot2::ggplot(x, aes(x = contraction, y = z_score)) + geom_point(alpha = alpha) +
    expand_limits(x = c(0,1)) +
    facet_wrap(~variable)
}


#' Plot the simulated "true" values versus posterior estimates
#'
#' @param x object containing results (a data.frame or [SBC_results()] object).
#' @param variables variables to show in the plot or `NULL` to show all
#' @param estimate which estimate to use for the central tendency,
#' must correspond a field already computed in the results (most likely `"mean"` and `"median"`).
#' @param uncertainty which estimates to use for uncertainty (a character vector of length 2)
#' must correspond a field already computed in the results. Pass `NULL` to avoid showing uncertainty at all.
#' @param alpha the alpha for the points and uncertainty intervals
#' @param parameters DEPRECATED, use `variables` instead
#' @return a ggplot2 plot object
#' @export
plot_sim_estimated <- function(x, variables = NULL, estimate = "mean",
                               uncertainty = c("q5", "q95"),
                               alpha = NULL, ..., parameters = NULL) {
  UseMethod("plot_sim_estimated")
}

#' @export
plot_sim_estimated.SBC_results <- function(x, variables = NULL, estimate = "mean",
                                           uncertainty = c("q5", "q95"),
                                           alpha = NULL, ..., parameters = NULL) {
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  plot_sim_estimated(x$stats, variables = variables, estimate = estimate,
                     uncertainty = uncertainty, alpha = alpha, ...)
}

#' @param show_hidden Show variables marked with [hidden_var_attribute()]
#'    (by default, those are not shown, available only when `x` is a data.frame)
#' @export
#' @rdname plot_sim_estimated
plot_sim_estimated.data.frame <- function(x, variables = NULL, estimate = "mean",
                                          uncertainty = c("q5", "q95"),
                                          alpha = NULL,
                                          show_hidden = FALSE,
                                          parameters = NULL) {

  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  # Ensuring backwards compatibility
  if("parameter" %in% names(x)) {
    if(!("variable" %in% names(x))) {
      warning("The x parameter contains a `parameter` column, which is deprecated, use `variable` instead.")
      x$variable <- x$parameter
    }
  }

  required_columns <- c("variable", estimate, uncertainty)
  if(!all(required_columns %in% names(x))) {
    stop("The data.frame needs to have the following columns: ", paste0("'", required_columns, "'", collapse = ", "))
  }

  x <- filter_stats_by_variables_and_hidden(x, variables, show_hidden)

  if(is.null(alpha)) {
    n_points <- dplyr::summarise(dplyr::group_by(x, variable), count = dplyr::n())
    max_points <- max(n_points$count)
    alpha_guess <- 1 / ((max_points * 0.06) + 1)
    alpha <-  max(0.05, alpha_guess)
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
    y_label <- paste0("posterior ", estimate, " (", uncertainty[1], " - ", uncertainty[2], ")")
  } else {
    main_geom <- geom_point(alpha = alpha)
    all_aes <- aes(x = simulated_value, y = estimate__)
    y_label <- paste0("posterior ", estimate)
  }

  if(nrow(x) == 0) {
    stop("No data to plot.")
  }

  ggplot2::ggplot(x, all_aes) +
    geom_abline(intercept = 0, slope = 1, color = "skyblue1", linewidth = 2) +
    main_geom +
    labs(y = y_label) +
    facet_wrap(~variable, scales = "free")
}


#' Plot the observed coverage and its uncertainty.
#'
#' `plot_coverage` will plot the observed coverage,
#' while `plot_coverage_diff` will show the difference between observed
#' and expected coverage.
#' Please refer to [empirical_coverage()] for details on computation
#' and limitations of this plot as well as details on the arguments.
#' See `vignette("rank_visualizations")` for
#' more details.
#'
#' @param x object containing results (a data.frame or [SBC_results()] object).
#' @param variables variables to show in the plot or `NULL` to show all
#' @param prob the with of the uncertainty interval to be shown
#' @param max_points maximum number of points where to evaluate the coverage.
#'   If set to `NULL`, coverage is evaluated across the whole range of ranks.
#'   Setting to some smaller number may reduce memory footprint and increase speed.
#' @param parameters DEPRECATED. Use `variables` instead.
#' @return a ggplot2 plot object
#' @seealso empirical_coverage
#' @export
plot_coverage <- function(x, variables = NULL, prob = 0.95,
                          interval_type = "central", ..., parameters = NULL,
                          max_points = NULL) {
  UseMethod("plot_coverage")
}

#' @export
plot_coverage.SBC_results <- function(x, variables = NULL, prob = 0.95,
                                      interval_type = "central", ...,
                                      parameters = NULL,
                                      max_points = NULL) {

  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  plot_coverage(x$stats, variables = variables, prob = prob, interval_type = interval_type,
                max_points = max_points, ...)
}

#' @export
plot_coverage.data.frame <- function(x, variables = NULL, prob = 0.95,
                                     interval_type = "central",
                                     show_hidden = FALSE,
                                     parameters = NULL,
                                     max_points = NULL) {

  # Ensuring backwards compatibility
  if("parameter" %in% names(x)) {
    if(!("variable" %in% names(x))) {
      warning("The x parameter contains a `parameter` column, which is deprecated, use `variable` instead.")
      x$variable <- x$parameter
    }
  }

  if(!all(c("variable", "rank", "max_rank") %in% names(x))) {
    stop(SBC_error("SBC_invalid_argument_error",
                   "The stats data.frame needs a 'variable', 'rank' and 'max_rank' columns"))
  }

  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }

  x <- filter_stats_by_variables_and_hidden(x, variables, show_hidden)

  max_max_rank <- max(x$max_rank)
  if(is.null(max_points) || max_max_rank + 2 <= max_points) {
    widths <- (0:(max_max_rank + 1)) / (max_max_rank + 1)
  } else {
    widths <- seq(0,1, length.out = max_points)
  }
  coverage <- empirical_coverage(x, widths, prob = prob,
                                interval_type = interval_type)

  ggplot2::ggplot(coverage, aes(x = width_represented, y = estimate,
                                ymin = ci_low, ymax = ci_high)) +
    geom_ribbon(fill = "black", alpha = 0.33) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, color = "skyblue1", linewidth = 2) +
    geom_line() +
    scale_x_continuous(paste0(interval_type, " interval width"),
                       labels = scales::percent) +
    scale_y_continuous("Observed coverage", labels = scales::percent) +
    facet_wrap(~variable)
}


#' @rdname plot_coverage
#' @export
plot_coverage_diff <- function(x, variables = NULL, prob = 0.95,
                          interval_type = "central",
                          ...,
                          parameters = NULL,
                          max_points = NULL) {
  UseMethod("plot_coverage_diff")
}

#' @export
plot_coverage_diff.SBC_results <- function(x, variables = NULL, prob = 0.95,
                                      interval_type = "central",
                                      ...,
                                      max_points = NULL) {
  plot_coverage_diff(x$stats, variables = variables, prob = prob,
                     interval_type = interval_type, max_points = max_points, ...)
}

#' @export
plot_coverage_diff.data.frame <- function(x, variables = NULL, prob = 0.95,
                                     interval_type = "central",
                                     show_hidden = FALSE,
                                     parameters = NULL,
                                     max_points = NULL) {
  if(!all(c("variable", "rank", "max_rank") %in% names(x))) {
    stop(SBC_error("SBC_invalid_argument_error",
                   "The stats data.frame needs a 'variable', 'rank' and 'max_rank' columns"))
  }

  x <- filter_stats_by_variables_and_hidden(x, variables, show_hidden)

  max_max_rank <- max(x$max_rank)
  if(is.null(max_points) || max_max_rank + 2 <= max_points) {
    widths <- (0:(max_max_rank + 1)) / (max_max_rank + 1)
  } else {
    widths <- seq(0,1, length.out = max_points)
  }
  coverage <- empirical_coverage(x, widths, prob = prob,
                                 interval_type = interval_type)

  coverage <- dplyr::mutate(coverage,
                            diff = estimate - width_represented,
                            diff_low = ci_low - width_represented,
                            diff_high = ci_high - width_represented
                            )

  ggplot2::ggplot(coverage, aes(x = width_represented, y = diff,
                                ymin = diff_low, ymax = diff_high)) +
    geom_ribbon(fill = "black", alpha = 0.33) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 0, color = "skyblue1", linewidth = 2) +
    geom_line() +
    scale_x_continuous(paste0(interval_type, " interval width"),
                       labels = scales::percent) +
    scale_y_continuous("Coverage diff", labels = scales::percent) +
    facet_wrap(~variable, scales = "free_y")
}
