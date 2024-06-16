#' @export
#' @rdname variable-attributes
var_attributes <- function(...) {
  var_attr <- structure(list(...), class = "var_attributes")
  return(validate_var_attributes(var_attr))
}

var_attributes_from_list <- function(attr_names, attr_list) {
  if(length(attr_names) == 0) {
    return(NULL)
  }
  stopifnot(length(attr_list) == 1 || length(attr_list) == length(attr_names))
  if(length(attr_list) == 1 && length(attr_names) > 1) {
    attr_list <- rep(attr_list, length(attr_names))
  }

  names(attr_list) <- attr_names
  do.call(var_attributes, attr_list)
}

#' @export
#' @rdname variable-attributes
validate_var_attributes <- function(var_attr) {
  if(is.null(var_attr)) {
    return(NULL)
  }

  if(!is.list(var_attr) || !inherits(var_attr, "var_attributes")) {
    stop("`var_attributes` must be a list")
  }

  if(is.null(names(var_attr)) || any(names(var_attr) == "")) {
    stop("All elements of var_attributes need to have non-empty names")
  }

  if(any(grepl("\\[|\\]", names(var_attr)))) {
    stop("Names of var_attributes must not contain brackets (attributes are assumed identical for all array elements)")
  }

  character_attributes <- purrr::map_lgl(var_attr, is.character)
  if(!all(character_attributes)) {
    message(paste0(names(var_attr)[!character_attributes], collapse = ", "))
    stop("All elements of `var_attributes` must be character vectors")
  }

  return(var_attr)
}

variable_names_to_var_attributes_names <- function(variable_names) {
  gsub(r"(\[.*\])", "", variable_names)
}

#' Predefined constants to use for variable attributes recognized by SBC.
#'
#' Attributes give additional information useful for presenting SBC results
#' concerning the variables.
#'
#' Should be passed via an instance of the `var_attributes` class to
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
#'
#' `na_valid_var_attribute` will treat NAs as potentially equal to any
#' other value in rank ordering. This gives the expected results when NA
#' represents rare problems in computation that
#' should be ignored (see [calculate_ranks_draws_matrix()]). Setting this
#' attribute also changes the ESS/Rhat computation to ignore NAs.
#'
#' `inf_valid_var_attribute` means infinity values may appear in the samples
#' (this is useful e.g. to note that the parameter is actually not present for the
#' given draw). Setting this
#' attribute also changes the ESS/Rhat computation to ignore infinities.
#'
#' `submodel_var_attribute` signals that the parameter belongs to a submodel
#' which can be extracted individually
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
na_valid_var_attribute <- function() {
  "na_valid"
}


#' @rdname variable-attributes
#' @export
inf_valid_var_attribute <- function() {
  "inf_valid"
}

#' @rdname variable-attributes
#' @export
submodel_var_attribute <- function(sub_id) {
  paste0("submodel(",sub_id,")")
}


#' Get attribute presence for a vector of variable names.
#'
#' @param attribute the attribute to test
#' @param variable_names names of variable to test for the attribute
#' @param var_attr an object defining the attributes, see [variable-attributes].
#' @return a logical vector the same length as `variable_names`
#' @export
attribute_present <- function(attribute, variable_names, var_attr) {
  var_attr_names <- variable_names_to_var_attributes_names(variable_names)

  var_attr_present <- purrr::map_lgl(var_attr, \(x) attribute %in% x)
  result <- rep(FALSE, length(variable_names))
  names(result) <- variable_names

  present_indices <- var_attr_names %in% names(var_attr)
  result[present_indices] <- var_attr_present[var_attr_names[present_indices]]
  return(result)
}

#' Check an attribute  in the `$stats` element of the results
#' for presence of an attribute.
#' @export
#'
attribute_present_stats <- function(attribute, x) {
  stopifnot(is.character(attribute))
  stopifnot(length(attribute) == 1)
  grepl(paste0("(^|,) *\\Q",attribute,"\\E *($|,)"), x)
}


var_attributes_to_attributes_column <- function(var_attr, variables) {
  var_attr_names <- variable_names_to_var_attributes_names(variables)
  missing_names <- setdiff(unique(var_attr_names),
                           names(var_attr))
  for(n in missing_names){
    var_attr[[n]] <- ""
  }
  attributes_collapsed <- purrr::map_chr(var_attr, \(x) paste0(x, collapse = ", "))
  names(attributes_collapsed) <- names(var_attr)
  res <- attributes_collapsed[var_attr_names]
  names(res) <- NULL
  return(res)
}


#' @details
#' It is currently by design that multiple copies of an attribute are kept
#'
#' @export
#' @param ... the individual [var_attributes()] objects to combine.
combine_var_attributes <- function(...) {
  all_var_attr <- purrr::map(list(...), validate_var_attributes)

  all_names_duplicates <- do.call(c, purrr::map(all_var_attr, names))
  all_names <- sort(unique(all_names_duplicates))

  collect_single_name <- function(name) {
    do.call(c, purrr::map(all_var_attr, \(x) x[[name]]))
  }

  res <- purrr::map(all_names, collect_single_name)
  names(res) <- all_names
  validate_var_attributes(structure(res, class = "var_attributes"))
}

#' Remove an attribute from a `data.frame` column with attributes.
#'
#' @details
#' By design, only one copy of an attribute is removed, if there are multiple copies
#'
#' @export
remove_attribute_from_stats <- function(attribute, x) {
  res <- trimws(sub(paste0("(^|,) *\\Q",attribute,"\\E *($|,)"), "\\1", x))
  res_no_trail <- trimws(sub(",$", "", res))
  return(res_no_trail)
}
