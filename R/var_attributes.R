#' Predefined constants to use for variable attributes recognized by SBC.
#'
#' Attributes give additional information useful for presenting SBC results
#' concerning the variables.
#'
#' Should be passed via a named list to the `var_attributes` argument of
#' [SBC_datasets()] (e.g. by returning a `$var_attributes` element from a
#' function passed to [SBC_generator_function()]).
#'
#' `possibly_constant_var_attribute` attribute signals that having all
#' posterior draws identical is possible and thus no warnings should be
#' made for the resulting NAs in rhat and ESS checks.
#'
#' `binary_var_attribute` marks the attribute as a binary variable (0 or 1)
#' and thus eligible for some special visualisations (TODO).
#'
#' `hidden_var_attribute` will hide the variable in default visualisations,
#' unless the variable is explicitly mentioned.
#'
#' `na_lowest_var_attribute` will treat NAs in the variable as the lowest
#' in rank ordering. This gives the expected results when NA is a special value
#' indicating e.g. that the variable is
#' not present at all in a given fit. (see [calculate_ranks_draws_matrix()])
#'
#' `na_valid_var_attribute` will treat NAs as potentially equal to any
#' other value in rank ordering. This gives the expected results when NA
#' represents rare problems in computation that
#' should be ignored (see [calculate_ranks_draws_matrix()]). Setting this
#' attribute also silences warning for NAs in ranks.
#'
#'
#' In SBC results, the attributes of a variable are summarised in the
#' `attributes` column of the `$stats` data.frame. Use [attribute_present_stats()]
#' to check for presence of an attribute there.
#'
#' @rdname variable-attributes
#' @export
possibly_constant_var_attribute <- function() {
  "possibly_constant"
}

#' @rdname variable-attributes
#' @export
binary_var_attribute <- function() {
  "binary"
}

#' @rdname variable-attributes
#' @export
hidden_var_attribute <- function() {
  "hidden"
}

#' @rdname variable-attributes
#' @export
na_lowest_var_attribute <- function() {
  "na_lowest"
}

#' @rdname variable-attributes
#' @export
na_valid_var_attribute <- function() {
  "na_valid"
}

#' Get attribute presence for a vector of variable names.
#'
#' @return a logical vector the same length as `variable_names`
#' @export
attribute_present <- function(attribute, variable_names, var_attributes) {
  var_attributes_present <- purrr::map_lgl(var_attributes, \(x) attribute %in% x)
  result <- rep(FALSE, length(variable_names))
  names(result) <- variable_names

  present_indices <- variable_names %in% names(var_attributes)
  result[present_indices] <- var_attributes_present[names(result)[present_indices]]
  return(result)
}

#' Check an attribute  in the `$stats` element of the results
#' for presence of an attribute.
#' @export
#'
attribute_present_stats <- function(attribute, x) {
  grepl(paste0("(^|,) *",attribute," *($|,)"), x)
}


var_attributes_to_attributes_column <- function(var_attributes, variables) {
  missing_names <- setdiff(variables, names(var_attributes))
  for(n in missing_names){
    var_attributes[[n]] <- ""
  }
  attributes_collapsed <- purrr::map_chr(var_attributes, \(x) paste0(x, collapse = ", "))
  names(attributes_collapsed) <- names(var_attributes)
  res <- attributes_collapsed[variables]
  names(res) <- NULL
  return(res)
}
