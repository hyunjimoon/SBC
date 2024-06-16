# Note: needs to be kept in sync with combine_var_attributes_for_bf
combine_draws_matrix_for_bf <- function(dm0, dm1, model_draws, NA_raw_dm = FALSE, model_var = "model") {
  stopifnot(all(model_draws %in% c(0,1)))
  stopifnot(posterior::is_draws_matrix(dm0))
  stopifnot(posterior::is_draws_matrix(dm1))
  stopifnot(is.character(model_var) & length(model_var) == 1)

  stopifnot(posterior::ndraws(dm0) == posterior::ndraws(dm1))
  stopifnot(posterior::ndraws(dm0) == length(model_draws))

  dm0_raw <- dm0
  if(NA_raw_dm) {
    dm0_raw[model_draws == 1,] <- NA
  }
  posterior::variables(dm0_raw) <- paste0(".m0.",posterior::variables(dm0))

  dm1_raw <- dm1
  if(NA_raw_dm) {
    dm1_raw[model_draws == 0,] <- NA
  }
  posterior::variables(dm1_raw) <- paste0(".m1.",posterior::variables(dm1))


  dm_model <- posterior::draws_matrix(model = model_draws)
  posterior::variables(dm_model) <- model_var

  all_variables_names <- unique(c(posterior::variables(dm0), posterior::variables(dm1)))
  list_all_vars <- list()

  param_not_present <- -Inf

  for(v in all_variables_names) {
    if(v %in% posterior::variables(dm0)) {
      var0 <- as.numeric(dm0[, v])
    } else {
      var0 <- param_not_present
    }
    if(v %in% posterior::variables(dm1)) {
      var1 <- as.numeric(dm1[, v])
    } else {
      var1 <- param_not_present
    }
    list_all_vars[[v]] <- dplyr::if_else(model_draws == 0, var0, var1)
  }
  all_vars <- do.call(posterior::draws_matrix, list_all_vars)
  posterior::bind_draws(
    dm_model,
    all_vars,
    dm0_raw,
    dm1_raw)
}

# Note: needs to be kept in sync with combine_draws_matrix_for_bf
combine_var_attributes_for_bf <- function(dm0, dm1, var_attr0, var_attr1, model_var = "model") {
  stopifnot(posterior::is_draws_matrix(dm0))
  stopifnot(posterior::is_draws_matrix(dm1))

  var_attr0 <- validate_var_attributes(var_attr0)
  var_attr1 <- validate_var_attributes(var_attr1)

  stopifnot(is.character(model_var) & length(model_var) == 1)


  raw_attrs <- function(dm, orig_attr, model_id) {
    attr_names <- variable_names_to_var_attributes_names(posterior::variables(dm))

    new_attr_vec <- c(hidden_var_attribute(), submodel_var_attribute(model_id))
    new_attr <- var_attributes_from_list(attr_names, list(new_attr_vec))

    attr_combined <- combine_var_attributes(new_attr, orig_attr)
    prefix <- paste0(".m", model_id, ".")
    names(attr_combined) <- paste0(prefix, names(attr_combined))

    return(attr_combined)
  }

  raw_attr0 <- raw_attrs(dm0, var_attr0, 0)
  raw_attr1 <- raw_attrs(dm1, var_attr1, 1)

  single_model_attr_names <- variable_names_to_var_attributes_names(c(
    setdiff(posterior::variables(dm0), posterior::variables(dm1)),
    setdiff(posterior::variables(dm1), posterior::variables(dm0))
  ))
  single_model_attr <- var_attributes_from_list(single_model_attr_names, list(inf_valid_var_attribute()))

  attr_model <- var_attributes(model = c(binary_var_attribute(), possibly_constant_var_attribute()))
  names(attr_model) <- model_var

  return(
    combine_var_attributes(attr_model,
                           single_model_attr,
                           var_attr0,
                           var_attr1,
                           raw_attr0,
                           raw_attr1
    )
  )
}

#' @export
SBC_datasets_for_bf <- function(datasets_H0, datasets_H1, prob_H1 = 0.5, model_var = "model") {
  validate_SBC_datasets(datasets_H0)
  validate_SBC_datasets(datasets_H1)
  stopifnot(length(datasets_H0) == length(datasets_H1))

  model_draws <- rbinom(length(datasets_H0), size = 1, prob = prob_H1)
  combined_variables <- combine_draws_matrix_for_bf(datasets_H0$variables, datasets_H1$variables,
                                                    model_draws,
                                                    NA_raw_dm = TRUE,
                                                    model_var = model_var)

  combined_var_attributes <- combine_var_attributes_for_bf(datasets_H0$variables, datasets_H1$variables,
                                                           datasets_H0$var_attributes, datasets_H1$var_attributes,
                                                           model_var = model_var)

  combined_generated <- datasets_H0$generated
  for(i in 1:length(datasets_H0)) {
    if(model_draws[i] == 1) {
      combined_generated[[i]] <- datasets_H1$generated[[i]]
    }
  }

  SBC_datasets(combined_variables, combined_generated, combined_var_attributes)
}


#' Extract stats for a submodel from a Bayes factor results
#' @export
get_stats_for_submodel <- function(stats, submodel_id) {
  stats_subm <- dplyr::filter(stats, attribute_present_stats(submodel_var_attribute(submodel_id), attributes))
  if(nrow(stats_subm) == 0) {
    warning(paste0("There were no variables marked with `submodel_var_attribute(", submodel_id,")`"))
  }
  prefix_regex <- paste0("^\\.m", submodel_id, "\\.")
  stats_subm <- dplyr::mutate(stats_subm,
                            variable = gsub(prefix_regex, "", variable),
                            attributes = remove_attribute_from_stats(hidden_var_attribute(), attributes),
                            attributes = remove_attribute_from_stats(submodel_var_attribute(submodel_id), attributes)
  )
  return(stats_subm)
}


#' Split stats from SBC of Bayes factor into stats for the individual submodels.
#'
#' @return list of two stats `data.frame`s, each for simulations where the
#' corresponding model was chosen.
#' @export
split_SBC_results_for_bf <- function(results_bf) {
  stats <- results_bf$stats
  list(
    stats_H0 = get_stats_for_submodel(stats, 0),
    stats_H1 = get_stats_for_submodel(stats, 1)
  )
}
