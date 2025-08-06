#' @export
binary_probabilities_from_stats <- function(stats) {
  stats <- dplyr::filter(stats, attribute_present_stats(binary_var_attribute(), attributes))
  if(nrow(stats) == 0) {
    warning("No variables annotated with `binary_var_attribute()` found.")
    return(data.frame(prob = c(), estimate = c(), low = c(), high()))
  }
  if(!(("cdf_low" %in% names(stats)) || !("cdf_high" %in% names(stats)))) {
    stop("The range of continous ranks (cdf_low and cdf_high) column must be present in stats.\nIf it is not present, it likely means your backend
         does not implement `SBC_posterior_cdf()`.")
  }
  if(!all(stats$cdf_low == 0 | stats$cdf_high == 1)) {
    stop("For binary variables either the cdf_low needs to be 0 or the cdf_high needs to be 1")
  }

  stats <- dplyr::mutate(stats, prob = if_else(cdf_low == 0, cdf_high, cdf_low))

  return(stats)
}

#' @export
binary_calibration_from_stats <- function(stats, type = c("reliabilitydiag", "calibrationband"), ...) {
  stats <- binary_probabilities_from_stats(stats)

  stats_grouped <- dplyr::group_by(stats, variable)
  res <- dplyr::reframe(stats_grouped, binary_calibration_base(prob, simulated_value, type = type, ...))

  return(res)
}

#' @export
binary_calibration_base <- function(prob, outcome, uncertainty_prob = 0.95, type = c("reliabilitydiag", "calibrationband"), ...) {
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
    require_package_version("reliabilitydiag", "0.2.1", "to compute binary calibration with the type 'reliabilitydiag'.")
    rel_diag <- reliabilitydiag::reliabilitydiag(
      x = prob,
      y = outcome,
      region.level = uncertainty_prob,
      ...
    )
    res <- data.frame(prob = rel_diag$x$regions$x, low = rel_diag$x$regions$lower, high = rel_diag$x$regions$upper)
    res$estimate <- approx(x = c(rel_diag$x$bins$x_min, rel_diag$x$bins$x_max),
                           y = rep(rel_diag$x$bins$CEP_pav, times = 2),
                           xout = res$prob)$y

    return(res)
  } else if(type == "calibrationband") {
    require_package_version("calibrationband", "0.2", "to compute binary calibration with the type 'calibrationband'.")
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

    return(res)
  } else {
    stop("Unknown type")
  }
}


#' @export
plot_binary_calibration_diff <- function(stats, type = c("reliabilitydiag", "calibrationband"), ...) {
  calib_df <- binary_calibration_from_stats(stats, type = type, ...)

  ggplot(calib_df, aes(x = prob, ymin = low - prob, ymax = high - prob, y = estimate - prob)) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 0, color = "skyblue1", size = 2) +
    geom_ribbon(fill = "black", alpha = 0.33) +
    geom_line() + facet_wrap(~variable)
}

#' @export
plot_binary_calibration <- function(stats, type = c("reliabilitydiag", "calibrationband"), ...) {
  calib_df <- binary_calibration_from_stats(stats, type = type, ...)

  ggplot(calib_df, aes(x = prob, ymin = low, ymax = high, y = estimate)) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, color = "skyblue1", size = 2) +
    geom_ribbon(fill = "black", alpha = 0.33) +
    geom_line() + facet_wrap(~variable) + coord_fixed()
}
