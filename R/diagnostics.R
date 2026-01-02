#' @export
get_diagnostic_messages_single <- function(diagnostic, values) {
  UseMethod("get_diagnostic_messages_single")
}

#' @export
numeric_diagnostic <- function(label, report = NULL, error_above = Inf, error_below = -Inf, allow_na = FALSE, label_short = NULL,
                                   hint = "", digits = 2, unit = "") {
  if(is.null(label_short)) {
    label_short <- label
  }

  structure(list(label = label, report = report, error_above = error_above,
                 error_below = error_below, allow_na = allow_na,
                 label_short = label_short, hint = hint, digits = digits,
                 unit = unit),
            class = "SBC_numeric_diagnostic")
}


#' @export
get_diagnostic_messages_single.SBC_numeric_diagnostic <- function(diagnostic, values) {
  stopifnot(is.numeric(values) || is.integer(values))

  unit <- diagnostic$unit
  if(is.finite(diagnostic$error_below)) {
    if(is.finite(diagnostic$error_above)) {
      stop("Not implemented yet")
    } else {
      threshold_summary_error <- paste0("< ", diagnostic$error_below, unit)
      threshold_summary_ok <- paste0(">= ", diagnostic$error_below, unit)
    }
  } else {
    if(is.finite(diagnostic$error_above)) {
      threshold_summary_error <- paste0("> ", diagnostic$error_above, unit)
      threshold_summary_ok <- paste0("<= ", diagnostic$error_above, unit)
    } else {
      threshold_summary_error <- NULL
      threshold_summary_ok <- NULL
    }
  }

  if(is.null(threshold_summary_ok)) {
    label_for_report <- diagnostic$label
  } else {
    label_for_report <- diagnostic$label_short
  }

  if(is.null(diagnostic$report)) {
    report <- ""
  } else if(length(values) == 0) {
    report <- "No values available"
  } else if(diagnostic$report == "max") {
    report <- paste0("Maximum ", label_for_report, " was ", round(max(values, na.rm = TRUE), digits = diagnostic$digits), unit, ". ")
  } else if(diagnostic$report == "min") {
    report <- paste0("Minimum ", label_for_report, " was ", round(min(values, na.rm = TRUE), digits = diagnostic$digits), unit, ". ")
  } else if(diagnostic$report == "quantiles") {
    quantiles <- paste0(round(quantile(values, probs = c(0.05, 0.5, 0.95), na.rm = TRUE), digits = diagnostic$digits), unit)
    report <- paste0(label_for_report, " Q.05 = ", quantiles[1],", median = ", quantiles[2], ", Q.95 = ", quantiles[3], ". ")
  } else {
    stop(paste0("Unrecognized 'report' value: ", report))
  }


  ok <-
    !any(values > diagnostic$error_above, na.rm = TRUE) &&
    !any(values < diagnostic$error_below, na.rm = TRUE)

  n_na <- sum(is.na(values))

  if(is.null(threshold_summary_ok)) {
    if(!diagnostic$allow_na && n_na > 0) {
      msg <- data.frame(type = "bad", message = paste0(n_na, " (", round(100 * n_na / length(values)), "%) fits had NA ", label, report))
    } else {
      msg <- data.frame(type = "info", message = report)
    }
  } else {
    if(!diagnostic$allow_na && n_na > 0) {
      ok <- FALSE
      threshold_summary_error <- paste0(threshold_summary_error, " or NA")
    }

    if(ok) {
      msg <- data.frame(type = "ok", message = paste0("All fits had ", diagnostic$label, " ", threshold_summary_ok, ". ", report))
    } else {
      n_wrong <- sum(values > diagnostic$error_above | values < diagnostic$error_below, na.rm = TRUE)
      msg <- data.frame(type = "bad",
                        message = paste0( n_wrong, " (", round(100 * n_wrong / length(values)), "%) fits had ",
                                          diagnostic$label, " ", threshold_summary_error, ". ", report, diagnostic$hint))
    }
  }
  msg
}


#' @export
count_diagnostic <- function(label, error_above = 0, label_short = NULL, hint = "",
                             error_only = FALSE) {
  if(is.null(label_short)) {
    label_short <- label
  }
  structure(list(label = label,
                 error_above = error_above,
                 label_short = label_short,
                 hint = hint,
                 error_only = error_only),
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
    ok_type <- "info"
    threshold_summary <- "some "
  } else {
    ok_type <- "ok"
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
    msg <- data.frame(type = "bad",
                      message = paste0( n_above, " (", round(100 * n_above / length(values)), "%) fits had ",
                                        threshold_summary, diagnostic$label, ". ", report, diagnostic$hint))
  } else if(diagnostic$error_only) {
    msg <- data.frame(type = character(),
                      message = character())
  } else {
    msg <- data.frame(type = ok_type,
                      message = ok_message)
  }

  msg
}

#' @export
logical_diagnostic <- function(ok_value, true_label,
                               false_label = paste0("not ", true_label),
                               hint = "", error_only = FALSE) {
  stopifnot(!is.na(ok_value))
  structure(list(true_label = true_label,
                 false_label = false_label,
                 ok_value = ok_value,
                 error_only = error_only,
                 hint = hint),
            class = "SBC_logical_diagnostic")
}

#' @export
get_diagnostic_messages_single.SBC_logical_diagnostic <- function(diagnostic, values) {
  stopifnot(is.logical(values))
  ok_vec_no_nas <- tidyr::replace_na(values, !diagnostic$ok_value) == diagnostic$ok_value

  if(all(ok_vec_no_nas)) {
    if(diagnostic$error_only) {
       data.frame(type = character(), message = character())
    } else if(diagnostic$ok_value) {
      data.frame(type = "ok", message = paste0("All fits ", diagnostic$true_label,"."))
    } else {
      data.frame(type = "ok", message = paste0("No fits ", diagnostic$true_label,"."))
    }
  } else {
    n_fail <- sum(!ok_vec_no_nas)
    if(diagnostic$ok_value) {
      fail_label <- diagnostic$false_label
    } else {
      fail_label <- diagnostic$true_label
    }
    data.frame(type = "bad",
                      message = paste0( n_fail, " (", round(100 * n_fail / length(values)), "%) fits ",
                                        fail_label,". ", diagnostic$hint))
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
    data.frame(type = NA, message =
                 paste0(diagnostic$col_name, ": ", n_true, " (", round(100 * n_true / length(values)), "%) true, ",
                        n_na, " (", round(100 * n_na / length(values)), "%) NA, "))
  } else if(is.numeric(values) || is.integer(values)) {
    n_na <- sum(is.na(values))
    data.frame(type = NA, message =
                 paste0(diagnostic$col_name, ": min = ", round(min(values, na.rm = TRUE), digits = 3), ", max = ", round(max(values, na.rm = TRUE, digits = 3)), ", ",
                        n_na, " (", round(100 * n_na / length(values)), "%) NA, "))

  } else{
    tab <- sort(table(values, useNA = "ifany"), decreasing = TRUE)
    if(length(tab) <= 4) {
      tab_str <- paste0(names(tab), " - ", tab, " (", round(100 * tab / length(values)), "%)", collapse = ", ")
    } else {
      tab_str <- paste0(length(tab), " values. Most common: ", names(tab)[1], " - ", tab[1], " (", round(100 * tab[1] / length(values)), "%)")
    }
    data.frame(type = NA,
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
  data.frame(type = character(), message = character())
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
  if(!all(is.na(msgs$type) | (msgs$type %in% c("info", "ok", "bad")))) {
    warning("Urecognized message types")
  }
  msgs
}
