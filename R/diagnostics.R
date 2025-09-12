#' @export
SBC_numeric_diagnostic <- function(label, report = NULL, lower_threshold = Inf, upper_threshold = -Inf, allow_na = FALSE, label_short = NULL,
                                   hint = "") {
  if(is.null(label_short)) {
    label_short <- label
  }

  structure(list(label = label, report = report, lower_threshold = lower_threshold,
                 upper_threshold = upper_threshold,
                 label_short = label_short, hint = hint),
            class = "SBC_numeric_diagnostic")
}

#' @export
SBC_count_diagnostic <- function(label, lower_threshold = 0, label_short = NULL, hint = "") {
  if(is.null(label_short)) {
    label_short <- label
  }
  structure(list(label = label,
                 lower_threshold = lower_threshold,
                 label_short = label_short,
                 hint = hint),
            class = "SBC_count_diagnostic")
}


#' @export
SBC_logical_diagnostic <- function(label, error_value, hint) {
   stop("Baf")
}

#' @export
SBC_submodel_diagnostic <- function(prefix, diag) {
  stop("Baf")
}

#' @export
SBC_get_diagnostic_messages <- function(diagnostic, values) {
  UseMethod("SBC_get_diagnostic_messages")
}



#' @export
SBC_get_diagnostic_messages.SBC_numeric_diagnostic <- function(diagnostic, values) {
  stopifnot(is.numeric(values) || is.integer(values))


  if(is.finite(diagnostic$lower_threshold)) {
    if(is.finite(diagnostic$upper_threshold)) {
      stop("Not implemented yet")
    } else {
      threshold_summary = paste0("> ", diagnostic$lower_threshold)
    }
  } else {
    if(is.finite(diagnostic$upper_threshold)) {
      threshold_summary = paste0("< ", diagnostic$upper_threshold)
    } else {
      threshold_summary = NULL
    }
  }

  if(is.null(threshold_summary)) {
    label_for_report <- diagnostic$label
  } else {
    label_for_report <- diagnostic$label_short
  }

  if(is.null(diagnostic$report)) {
    report <- ""
  } else if(diagnostic$report == "max") {
    report <- paste0("Maximum ", label_for_report, " was ", max(values), ". ")
  } else if(diagnostic$report == "min") {
    report <- paste0("Minimum ", label_for_report, " was ", min(values), ". ")
  }


  ok <-
    !any(values > diagnostic$lower_threshold, na.rm = TRUE) &&
    !any(values < diagnostic$upper_threshold, na.rm = TRUE)

  if(is.null(threshold_summary)) {
    msg <- data.frame(ok = TRUE, message = report)
  } else {
    if(ok) {
      msg <- data.frame(ok = TRUE, paste0("All fits had ", diagnostic$label, " ", threshold_summary, ". ", report))
    } else {
      n_wrong <- sum(values > diagnostic$lower_threshold | values < diagnostic$upper_threshold, na.rm = TRUE)
      msg <- data.frame(ok = FALSE, paste0( n_wrong, " (", round(100 * n_wrong / length(values)), "%) fits had ", diagnostic$label, " ", threshold_summary, ". ", report, diagnostic$hint))
    }
  }
  msg
}

#' @export
SBC_get_diagnostic_messages.SBC_count_diagnostic <- function(diagnostic, values) {
  stopifnot(is.numeric(values) || is.integer(values))
  stopifnot(all(as.integer(values) == values))

  any_above_thresh <- any(values > diagnostic$lower_threshold, na.rm = TRUE)

  if(diagnostic$lower_threshold == 0) {
    threshold_summary <- "some"
  } else {
    threshold_summary <- paste0("> ", diagnostic$lower_threshold)
  }

  if(any_above_thresh) {
    n_above <- sum(values > diagnostic$lower_threshold, na.rm = TRUE)
    msg <- data.frame(ok = FALSE,
                      message = paste0( n_above, " (", round(100 * n_above / length(values)), "%) fits had ",
                                        threshold_summary," ", diagnostic$label, ". Maximum number of ", diagnostic$label_short, " was ", max(values, na.rm = TRUE), ". ", diagnostic$hint))
  } else {
    msg <- data.frame(ok = TRUE,
                      message = paste0( "No fits had ",
                                        threshold_summary," ", diagnostic$label, "."))
  }

  msg
}

#' @export
SBC_get_diagnostic_message.SBC_logical_diagnostic <- function(diagnostic, values) {
  stopifnot(is.logical(values))
  stop("BFDFD")
}


#' @export
SBC_get_diagnostic_message.SBC_submodel_diagnostic <- function(diagnostic, values) {
  stop("TODO")
}


