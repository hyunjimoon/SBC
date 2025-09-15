#' @export
get_diagnostic_messages_single <- function(diagnostic, values) {
  UseMethod("get_diagnostic_messages_single")
}

#' @export
numeric_diagnostic <- function(label, report = NULL, error_above = Inf, error_below = -Inf, allow_na = FALSE, label_short = NULL,
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
get_diagnostic_messages_single.SBC_numeric_diagnostic <- function(diagnostic, values) {
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
  } else if(diagnostic$report == "quantiles") {
    quantiles <- round(quantile(values, probs = c(0.05, 0.5, 0.95)), digits = diagnostic$digits)
    report <- paste0(label_for_report, " Q.05 = ", quantiles[1],", median = ", quantiles[2], ", Q.95 = ", quantiles[3], ". ")
  } else {
    stop(paste0("Unrecognized 'report' value: ", report))
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
count_diagnostic <- function(label, error_above = 0, label_short = NULL, hint = "") {
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
get_diagnostic_messages_single.SBC_count_diagnostic <- function(diagnostic, values) {
  stopifnot(is.numeric(values) || is.integer(values))
  stopifnot(all(as.integer(values) == values))

  any_above_thresh <- any(values > diagnostic$error_above, na.rm = TRUE)

  report <- paste0("Maximum number of ", diagnostic$label_short, " was ", max(values, na.rm = TRUE), ". ")

  if (diagnostic$error_above == Inf) {
    ok_message <- report
    threshold_summary <- "some "
  } else {
    if (diagnostic$error_above == 0) {
      threshold_summary <- ""
      ok_message <- paste0( "No fits had ",
              threshold_summary, diagnostic$label, ".")
    } else {
      threshold_summary <- paste0("> ", diagnostic$error_above, " ")
      ok_message <- paste0( "No fits had ",
              threshold_summary, diagnostic$label, ". ", report)
    }
  }


  if(any_above_thresh) {
    n_above <- sum(values > diagnostic$error_above, na.rm = TRUE)
    msg <- data.frame(ok = FALSE,
                      message = paste0( n_above, " (", round(100 * n_above / length(values)), "%) fits had ",
                                        threshold_summary, diagnostic$label, ". ", report, diagnostic$hint))
  } else {
    msg <- data.frame(ok = TRUE,
                      message = ok_message)
  }

  msg
}

#' @export
logical_diagnostic <- function(error_label, error_value, hint = "") {
  structure(list(error_label = error_label,
                 error_value = error_value,
                 hint = hint),
            class = "SBC_logical_diagnostic")
}

#' @export
get_diagnostic_messages_single.SBC_logical_diagnostic <- function(diagnostic, values) {
  stopifnot(is.logical(values))
  ok_vec <- (values != diagnostic$error_value)
  if(all(ok_vec)) {
    msg <- data.frame(ok = TRUE, message = paste0("No fits ", diagnostic$error_label,"."))
  } else {
    n_fail <- sum(!ok_vec)
    msg <- data.frame(ok = FALSE,
                      message = paste0( n_fail, " (", round(100 * n_fail / length(values)), "%) fits ",
                                        diagnostic$error_label,". ", diagnostic$hint))
  }
}

#' @export
submodel_diagnostic <- function(prefix, diag) {
  structure(list(prefix = prefix,
                 diag = diag),
            class = "SBC_submodel_diagnostic")
}

#' @export
get_diagnostic_messages_single.SBC_submodel_diagnostic <- function(diagnostic, values) {
  msgs <- get_diagnostic_messages_single(diagnostic$diag, values)
  msgs$message <- paste0(diagnostic$prefix,": ", msgs$message)
  msgs
}


#' @export
default_diagnostic <- function(col_name) {
  structure(list(col_name = col_name),
            class = "SBC_default_diagnostic")
}

#' @export
get_diagnostic_messages_single.SBC_default_diagnostic <- function(diagnostic, values) {
  if(is.logical(values)) {
    n_true <- sum(values, na.rm = TRUE)
    n_na <- sum(is.na(values))
    data.frame(ok = NA, message =
                 paste0(diagnostic$col_name, ": ", n_true, " (", round(100 * n_true / length(values)), "%) true, ",
                        n_na, " (", round(100 * n_na / length(values)), "%) NA, "))
  } else if(is.numeric(values) || is.integer(values)) {
    n_na <- sum(is.na(values))
    data.frame(ok = NA, message =
                 paste0(diagnostic$col_name, ": min = ", round(min(values, na.rm = TRUE), digits = 3), ", max = ", round(max(values, na.rm = TRUE, digits = 3)), ", ",
                        n_na, " (", round(100 * n_na / length(values)), "%) NA, "))

  } else{
    tab <- sort(table(values, useNA = "ifany"), decreasing = TRUE)
    if(length(tab) <= 4) {
      tab_str <- paste0(names(tab), " - ", tab, " (", round(100 * tab / length(values)), "%)", collapse = ", ")
    } else {
      tab_str <- paste0(length(tab), " values. Most common: ", names(tab)[1], " - ", tab[1], " (", round(100 * tab[1] / length(values)), "%)")
    }
    data.frame(ok = NA,
               message = paste0(diagnostic$col_name, ": ", tab_str))
  }
}

# A diagnostic that is ignored for summaries and messages
#' @export
skip_diagnostic <- function() {
  structure(list(), class = "SBC_skip_diagnostic")
}

#' @export
get_diagnostic_messages_single.SBC_skip_diagnostic <- function(diagnostic, values) {
  data.frame(ok = logical(), message = character())
}


#' @export
get_all_diagnostic_messages <- function(diags, types) {
  shared_names <- intersect(names(diags), names(types))

  missing_types <- names(diags)[!(names(diags) %in% c("sim_id", names(types)))]
  types_default <- purrr::map(missing_types, default_diagnostic)
  names(types_default) <- missing_types

  missing_diags <- names(types)[!(names(types) %in% names(diags))]
  if(length(missing_diags) > 0) {
    warning("Following declared diacnostic types have no values available: ", paste0(missing_diags, collapse = ", "))
  }

  diags_used <- c(diags[shared_names], diags[missing_types])
  types_used <- c(types[shared_names], types_default)

  msgs <- purrr::map2_dfr(types_used, diags_used, get_diagnostic_messages_single)
  msgs
}
