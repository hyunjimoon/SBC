#' @title Extract binary probabilities from SBC stats.
#'
#' @description Takes all variables marked with [binary_var_attribute()] in the
#' statistics (the `$stats` field of an `SBC_results` object)
#' and transforms the reported CDF values into probabilities
#' that the underlying binary variable is 1.
#'
#' @details
#' To support exact binary probabilities, a backend must implement [SBC_posterior_cdf()],
#' as exact probabilities cannot be derived from samples.
#'
#'
#' @returns A `data.frame` with at least the columns `variable` and `prob`.
#'
#' @export
binary_probabilities_from_stats <- function(stats) {
  stats <- dplyr::filter(stats, attribute_present_stats(binary_var_attribute(), attributes))
  if(nrow(stats) == 0) {
    warning("No variables annotated with `binary_var_attribute()` found.")
    return(data.frame(variable = c(), prob = c()))
  }
  if(!(("cdf_low" %in% names(stats)) || !("cdf_high" %in% names(stats)))) {
    stop("The range of continous ranks (cdf_low and cdf_high) column must be present in stats.\nIf it is not present, it likely means your backend
         does not implement `SBC_posterior_cdf()`.")
  }
  if(!("simulated_value" %in% names(stats))) {
    stop("The simulated_value column is required in stats")
  }
  if(!all(stats$cdf_low == 0 | stats$cdf_high == 1)) {
    stop("For binary variables either the cdf_low needs to be 0 or the cdf_high needs to be 1")
  }

  stats <- dplyr::mutate(stats, prob = dplyr::case_when(
    cdf_low == 0 & cdf_high == 1 ~ simulated_value,
    cdf_low == 0 ~ 1 - cdf_high,
    TRUE ~ 1 - cdf_low))

  return(stats)
}

#' Binary prediction calibration computation and visualisation.
#'
#' @description
#' Obtain estimate of binary prediction calibration and the associated
#' uncertainty interval for further processing or plot it directly.
#'
#' @param bp the binary probabilities --- typically obtained with
#' [binary_probabilities_from_stats()]. Can however be manually constructed,
#' it needs to be a `data.frame` with columns `variable`, `prob` and `simulated_value`.
#' @param type the type of calibration uncertainty bands to compute, see details.
#' @param alpha the level associated with the confidence intervals reports
#' @param region.position for `type ="reliabilitydiag"` we may choose whether
#' the uncertainty interval surrounds the estimate (`region.position = "estimate"`)
#' or the null distribution (`region.position = "diagonal"`)
#' @param ... additional arguments passed to
#'  [reliabilitydiag::reliabilitydiag()] or [calibrationband::calibration_bands()]
#'
#' @details
#' When `type = "reliabilitydiag"`, the intervals are based on
#' [reliabilitydiag::reliabilitydiag()] and depending on `region.position`
#' can be centered on the null distribution of perfect calibration or on the
#' estimated calibration.
#' When `type = "calibrationband"` the intervals
#' are around the estimated calibration using [calibrationband::calibration_bands()]
#' --- in our experience the `calibrationband`
#' method has less sensitivity to detect miscalibration, but they require
#' somewhat weaker assumptions.
#'
#' @returns `binary_calibration_from_bp` returns a `data.frame` with columns `variable`, `prob`, `estimate`, `low` and `high`,
#' for each variable, it contains an estimate + confidence interval across a range
#' of probabilities (in equal steps). The `plot_` methods return a `ggplot2`
#' object showing either the calibration curve or the difference between
#' the calibration curve and perfect calibration (the diagonal)
#'
#'
#' @rdname binary_calibration
#' @export
binary_calibration_from_bp <- function(bp, type = c("reliabilitydiag", "calibrationband"), alpha = 0.05, ..., region.position = NULL) {

  bp_grouped <- dplyr::group_by(bp, variable)
  res <- dplyr::reframe(bp_grouped, binary_calibration_base(prob, simulated_value, type = type, uncertainty_prob = 1 - alpha, region.position = region.position, ...))

  max_sims <- max(dplyr::tally(bp_grouped)$n)
  attr(res, "bp") <- bp
  attr(res, "bins") <- max(2, min(100, ceiling(max_sims / 10)))
  return(res)
}

binary_calibration_base <- function(prob, outcome, uncertainty_prob = 0.95, type = c("reliabilitydiag", "calibrationband"), ..., region.position = NULL) {
  stopifnot(is.numeric(prob))
  stopifnot((is.numeric(outcome) || is.logical(outcome) || is.integer(outcome)))
  outcome <- as.numeric(outcome)
  stopifnot(all(outcome %in% c(0,1)))
  stopifnot(all(prob >=0 & prob <= 1))
  stopifnot(length(prob) == length(outcome))

  na_indices <- is.na(prob) | is.na(outcome)

  prob <- prob[!na_indices]
  outcome <- outcome[!na_indices]

  type <- match.arg(type)
  if(type == "reliabilitydiag") {
    list_args <- list(...)
    if(is.null(region.position)) {
      region.position <- "diagonal"
    }


    require_package_version("reliabilitydiag", "0.2.1", "to compute binary calibration with the type 'reliabilitydiag'.")
    rel_diag <- reliabilitydiag::reliabilitydiag(
      x = prob,
      y = outcome,
      region.level = uncertainty_prob,
      region.position = region.position,
      ...
    )
    res <- data.frame(prob = rel_diag$x$regions$x, low = rel_diag$x$regions$lower, high = rel_diag$x$regions$upper)
    res$estimate <- approx(x = c(rel_diag$x$bins$x_min, rel_diag$x$bins$x_max),
                           y = rep(rel_diag$x$bins$CEP_pav, times = 2),
                           xout = res$prob)$y

    if(region.position == "diagonal") {
      res$interval_type = "null"
    } else if(region.position == "estimate") {
      res$interval_type = "estimate"
    } else {
      stop("unrecognized region.position")
    }

    return(res)
  } else if(type == "calibrationband") {
    require_package_version("calibrationband", "0.2", "to compute binary calibration with the type 'calibrationband'.")

    if(!is.null(region.position) && region.position != "estimate") {
      stop("calibrationband only supports region.position = 'estimate'")
    }
    # Need to remove extreme indices because they cause crashes in the package
    extreme_indices <- prob < 1e-10 | prob > 1 - 1e-10
    extreme_indices_mismatch <- extreme_indices & round(prob) != outcome
    if(any(extreme_indices_mismatch)) {
      print(which(extreme_indices_mismatch))
      warning("Some probabilities are very close 0/1 but are not matched by a model 0/1")
    }
    prob <- prob[!extreme_indices]
    outcome <- outcome[!extreme_indices]

    # Avoiding https://github.com/marius-cp/calibrationband/issues/1
    prob <- round(prob, digits = 7)

    bands <- calibrationband::calibration_bands(prob, outcome, alpha = 1 - uncertainty_prob, ...)

    res <- dplyr::transmute(bands$bands, prob = x, low = lwr, high = upr)

    data_estimate_raw <- tidyr::pivot_longer(
      bands$bins,
      cols = dplyr::all_of(c("x_min", "x_max")),
      values_to = "x")
    estimate_x <- c(0, data_estimate_raw$x, 1)
    estimate_isoy <- c(0, data_estimate_raw$isoy, 1)

    res$estimate <- approx(estimate_x, estimate_isoy, xout = res$prob, ties = "mean")$y
    res$interval_type = "estimate"

    return(res)
  } else {
    stop("Unknown type")
  }
}


calibration_prob_hist_geom <- function(calib_df) {
  geom_histogram(
    data = attr(calib_df, "bp"),
    aes(x = prob, y = after_stat(width * density)), # Making each subplot sum to 1, following https://github.com/donaldtmcknight/Proportional-histograms-and-density-plots-in-ggplot/blob/main/Tutorial_git.md#511-grouped-histogram-where-each-group-sums-to-1
    inherit.aes = FALSE,
    bins = attr(calib_df, "bins"),
    boundary = 0,
    fill = "transparent",
    color = "orangered"
  )
}

#' @param prob_histogram Whether a histogram of the observed probabilities should
#' be overlaid with the calibration curve.
#' @rdname binary_calibration
#' @export
plot_binary_calibration_diff <- function(x, type = c("reliabilitydiag", "calibrationband"), alpha = 0.05, ..., region.position = NULL, prob_histogram = TRUE) {
   UseMethod("plot_binary_calibration_diff", x)
}

#' @param res An [SBC_results] object
#' @rdname binary_calibration
#' @export
plot_binary_calibration_diff.SBC_results <- function(res, ...) {
  plot_binary_calibration_diff(
    binary_probabilities_from_stats(res$stats)
    , ...)
}

#' @rdname binary_calibration
#' @export
plot_binary_calibration_diff.data.frame <- function(bp, type = c("reliabilitydiag", "calibrationband"), alpha = 0.05, ..., region.position = NULL, prob_histogram = TRUE) {
  calib_df <- binary_calibration_from_bp(bp, type = type, alpha = alpha, region.position = region.position, ...)

  if(prob_histogram) {
    hist_geom <- calibration_prob_hist_geom(calib_df)
  } else {
    hist_geom <- NULL
  }

  ggplot(calib_df, aes(x = prob, ymin = low - prob, ymax = high - prob, y = estimate - prob)) +
    hist_geom +
    geom_segment(x = 0, y = 0, xend = 1, yend = 0, color = "skyblue1", size = 2) +
    geom_ribbon(aes(fill = interval_type), alpha = 0.33) +
    scale_fill_manual(values = c("null" = "skyblue1", "estimate" = "black"), guide = "none") +
    geom_line() + facet_wrap(~variable)
}


#' @param prob_histogram Whether a histogram of the observed probabilities should
#' be overlaid with the calibration curve.
#' @rdname binary_calibration
#' @export
plot_binary_calibration <- function(x, type = c("reliabilitydiag", "calibrationband"), alpha = 0.05, ..., region.position = NULL, prob_histogram = TRUE) {
  UseMethod("plot_binary_calibration", x)
}

#' @param res An [SBC_results] object
#' @rdname binary_calibration
#' @export
plot_binary_calibration.SBC_results <- function(res, ...) {
  plot_binary_calibration(
    binary_probabilities_from_stats(res$stats)
    , ...)
}

#' @rdname binary_calibration
#' @export
plot_binary_calibration.data.frame <- function(bp, type = c("reliabilitydiag", "calibrationband"), ..., region.position = NULL, prob_histogram = TRUE) {
  calib_df <- binary_calibration_from_bp(bp, type = type, ..., region.position = region.position)

  if(prob_histogram) {
    hist_geom <- calibration_prob_hist_geom(calib_df)
  } else {
    hist_geom <- NULL
  }

  #TODO also plot the histogram
  ggplot(calib_df, aes(x = prob, ymin = low, ymax = high, y = estimate)) +
    hist_geom +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, color = "skyblue1", size = 2) +
    geom_ribbon(aes(fill = interval_type), alpha = 0.33) +
    scale_fill_manual(values = c("null" = "skyblue1", "estimate" = "black"), guide = "none") +
    geom_line() + facet_wrap(~variable) + coord_fixed()
}
