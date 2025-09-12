#' @export
SBC_numeric_diagnostic <- function(label, report = NULL, error_above = Inf, error_below = -Inf, allow_na = FALSE, label_short = NULL,
                                   hint = "", digits = 2) {
  if(is.null(label_short)) {
    label_short <- label
  }

  structure(list(label = label, report = report, error_above = error_above,
                 error_below = error_below,
                 label_short = label_short, hint = hint, digits = digits),
            class = "SBC_numeric_diagnostic")
}

#' @export
SBC_count_diagnostic <- function(label, error_above = 0, label_short = NULL, hint = "") {
  if(is.null(label_short)) {
    label_short <- label
  }
  structure(list(label = label,
                 error_above = error_above,
                 label_short = label_short,
                 hint = hint),
            class = "SBC_count_diagnostic")
}


#' @export
SBC_logical_diagnostic <- function(label, error_value, hint) {
  structure(list(label = label,
                 error_value = error_value,
                 hint = hint),
            class = "SBC_logical_diagnostic")
}

#' @export
SBC_submodel_diagnostic <- function(prefix, diag) {
  structure(list(prefix = prefix,
                 diag = diag),
            class = "SBC_submodel_diagnostic")
}

#' @export
SBC_get_diagnostic_messages <- function(diagnostic, values) {
  UseMethod("SBC_get_diagnostic_messages")
}



#' @export
SBC_get_diagnostic_messages.SBC_numeric_diagnostic <- function(diagnostic, values) {
  stopifnot(is.numeric(values) || is.integer(values))


  if(is.finite(diagnostic$error_above)) {
    if(is.finite(diagnostic$error_below)) {
      stop("Not implemented yet")
    } else {
      threshold_summary = paste0("< ", diagnostic$error_above)
    }
  } else {
    if(is.finite(diagnostic$error_below)) {
      threshold_summary = paste0("> ", diagnostic$error_below)
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
    report <- paste0("Maximum ", label_for_report, " was ", round(max(values), digits = diagnostic$digits), ". ")
  } else if(diagnostic$report == "min") {
    report <- paste0("Minimum ", label_for_report, " was ", round(min(values), digits = diagnostic$digits), ". ")
  }


  ok <-
    !any(values > diagnostic$error_above, na.rm = TRUE) &&
    !any(values < diagnostic$error_below, na.rm = TRUE)

  if(is.null(threshold_summary)) {
    msg <- data.frame(ok = TRUE, message = report)
  } else {
    if(ok) {
      msg <- data.frame(ok = TRUE, message = paste0("All fits had ", diagnostic$label, " ", threshold_summary, ". ", report))
    } else {
      n_wrong <- sum(values > diagnostic$error_above | values < diagnostic$error_below, na.rm = TRUE)
      msg <- data.frame(ok = FALSE, message = paste0( n_wrong, " (", round(100 * n_wrong / length(values)), "%) fits had ", diagnostic$label, " ", threshold_summary, ". ", report, diagnostic$hint))
    }
  }
  msg
}

#' @export
SBC_get_diagnostic_messages.SBC_count_diagnostic <- function(diagnostic, values) {
  stopifnot(is.numeric(values) || is.integer(values))
  stopifnot(all(as.integer(values) == values))

  any_above_thresh <- any(values > diagnostic$error_above, na.rm = TRUE)

  if(diagnostic$error_above == 0) {
    threshold_summary <- "some"
  } else {
    threshold_summary <- paste0("> ", diagnostic$error_above)
  }

  if(any_above_thresh) {
    n_above <- sum(values > diagnostic$error_above, na.rm = TRUE)
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
SBC_get_diagnostic_messages.SBC_logical_diagnostic <- function(diagnostic, values) {
  stopifnot(is.logical(values))
  ok_vec <- (values == diagnostic$ok_value)
  if(all(ok_vec)) {
    msg <- data.frame(ok = TRUE, message = paste0("No fits ", diagnostic$label,"."))
  } else {
    n_fail <- sum(!ok_vec)
    msg <- data.frame(ok = FALSE,
                      message = paste0( n_fail, " (", round(100 * n_fail / length(values)), "%) fits ",
                                        diagnostic$label,". ", diagnostic$hint))

  }
}


#' @export
SBC_get_diagnostic_messages.SBC_submodel_diagnostic <- function(diagnostic, values) {
  stop("TODO")
}


