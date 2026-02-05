# Note: needs to be kept in sync with combine_var_attributes_for_bf
combine_draws_matrix_for_bf <- function(dm_list, model_draws, NA_raw_dm = FALSE, model_var = "model") {
  stopifnot(length(dm_list) >= 2)
  stopifnot(all(model_draws %in% seq(0, length(dm_list) - 1)))
  stopifnot(all(purrr::map_lgl(dm_list, posterior::is_draws_matrix)))
  stopifnot(is.character(model_var) & length(model_var) == 1)

  dm_n_draws <- purrr::map_int(dm_list, posterior::ndraws)
  if(length(unique(dm_n_draws)) > 1) {
    stop("Not all elements of dm_list have the same number of draws")
  }
  stopifnot(dm_n_draws[1] == length(model_draws))

  dm_list <- unname(dm_list)

  process_single_dm_raw <- function(dm, index_p1) {
    if(posterior::nvariables(dm) == 0) {
      return(dm)
    } else {
      dm_raw <- dm
      if(NA_raw_dm) {
        dm_raw[model_draws + 1 != index_p1,] <- NA
      }
      posterior::variables(dm_raw) <- paste0(".m", index_p1 - 1, ".",posterior::variables(dm))
      return(dm_raw)
    }
  }
  dm_raw_list <- purrr::imap(dm_list, process_single_dm_raw)

  pairs_index_by_model <- matrix(nrow = length(model_draws), ncol = 2)
  pairs_index_by_model[, 1] <- 1:length(model_draws)
  pairs_index_by_model[, 2] <- model_draws + 1


  dm_model <- posterior::draws_matrix(model = model_draws)
  posterior::variables(dm_model) <- model_var
  if(length(dm_list) > 2) {
    dm_is_model_raw <- matrix(0, nrow = length(model_draws), ncol = length(dm_list))
    dm_is_model_raw[pairs_index_by_model] <- 1
    colnames(dm_is_model_raw) <- paste0("is_", model_var, 0:(length(dm_list) - 1))
    dm_model <- posterior::bind_draws(dm_model, posterior::as_draws_matrix(dm_is_model_raw))
  }

  all_variables_names <- unique(do.call(c, purrr::map(dm_list, posterior::variables)))
  if(is.null(all_variables_names)) {
    all_vars <- posterior::as_draws_matrix(matrix(nrow = length(model_draws), ncol = 0))
  } else {
    list_all_vars <- list()

    param_not_present <- -Inf

    for(v in all_variables_names) {
      var_matrix <- matrix(nrow = length(model_draws), ncol = length(dm_list))
      for(m in 1:length(dm_list)) {
        if(v %in% posterior::variables(dm_list[[m]])) {
          var_matrix[, m] <- as.numeric(dm_list[[m]][, v])
        } else {
          var_matrix[, m] <- param_not_present
        }
      }
      list_all_vars[[v]] <- var_matrix[pairs_index_by_model]
    }
    all_vars <- do.call(posterior::draws_matrix, list_all_vars)
  }

  all_dm <- c(list(dm_model, all_vars), dm_raw_list)
  do.call(posterior::bind_draws, all_dm)
}

# Note: needs to be kept in sync with combine_draws_matrix_for_bf
combine_var_attributes_for_bf <- function(dm_list, var_attr_list, model_var = "model") {
  stopifnot(length(dm_list) == length(var_attr_list))
  stopifnot(length(dm_list) >= 2)
  stopifnot(all(purrr::map_lgl(dm_list, posterior::is_draws_matrix)))


  var_attr_list <- purrr::map(var_attr_list, validate_var_attributes)

  stopifnot(is.character(model_var) & length(model_var) == 1)


  raw_attrs_single <- function(dm, orig_attr, model_id) {
    attr_names <- variable_names_to_var_attributes_names(posterior::variables(dm))

    new_attr_vec <- c(hidden_var_attribute(), submodel_var_attribute(model_id))
    new_attr <- var_attributes_from_list(attr_names, list(new_attr_vec))

    attr_combined <- combine_var_attributes(new_attr, orig_attr)
    prefix <- paste0(".m", model_id, ".")
    names(attr_combined) <- paste0(prefix, names(attr_combined))

    return(attr_combined)
  }

  raw_attr_list <- purrr::pmap(list(dm = dm_list, orig_attr = var_attr_list, model_id = 0:(length(dm_list) - 1)),
                               raw_attrs_single)

  shared_attr_names <- intersect(posterior::variables(dm_list[[1]]), posterior::variables(dm_list[[2]]))
  if(length(dm_list) > 2) {
    for(i in 3:length(dm_list)) {
      shared_attr_names <- intersect(shared_attr_names, posterior::variables(dm_list[[i]]))
    }
  }

  not_shared_attr_single <- function(dm) {
    setdiff(posterior::variables(dm), shared_attr_names)
  }

  not_shared_attr_names <- variable_names_to_var_attributes_names(unique(do.call(
    c, purrr::map(dm_list, not_shared_attr_single)
  )))

  not_shared_attr <- var_attributes_from_list(not_shared_attr_names, list(inf_valid_var_attribute()))

  if(length(dm_list) == 2) {
    attr_model <- var_attributes(model = c(binary_var_attribute(), possibly_constant_var_attribute()))
    names(attr_model) <- model_var
  } else {
    attrs_list <- c(list(possibly_constant_var_attribute()),
                    rep(list(
                      c(binary_var_attribute(), possibly_constant_var_attribute())
                    ), times = length(dm_list)))
    attr_names <- c(model_var, paste0("is_", model_var, 0:(length(dm_list) - 1)))
    attr_model <- var_attributes_from_list(attr_names, attrs_list)
  }

  all_attr_list <- c(list(attr_model, not_shared_attr), var_attr_list, raw_attr_list)

  return(
    do.call(combine_var_attributes, all_attr_list)
  )
}

#' Combine two or more datasets to check Bayes factor computation.
#'
#' @description
#' Merge two or more datasets of the same length, each generated by a different model. The "true" model
#' will be chosen randomly for each dataset.
#'
#' @details The merged dataset will keep track of all parameters and later allow
#' splitting results with [`split_SBC_results_for_bf()`], to this, a number of
#' submodel-specific variables will be stored, but marked as
#' [`hidden_var_attribute()`] to not clutter default visualisations.
#' @export
SBC_datasets_for_bf <- function(..., probs = NULL, model_var = "model", prob_H1 = NULL) {
  all_datasets <- purrr::map(list(...), validate_SBC_datasets)

  if(!is.null(prob_H1)) {
    stopifnot(length(all_datasets) == 2)
    stopifnot(is.null(probs))
    probs <- c(1 - prob_H1, prob_H1)
  }

  if(is.null(probs)) {
    probs <- rep(1 / length(all_datasets), times = length(all_datasets))
  } else {
    stopifnot(is.numeric(probs))
    stopifnot(length(probs) == length(all_datasets))
  }

  all_lengths <- purrr::map_dbl(all_datasets, length)
  n_sims <- unique(all_lengths)
  if(length(n_sims) != 1) {
    stop("Datasets have different sizes")
  }

  if(length(all_datasets) == 2) {
    # Keep old code for this case to prevent invalidation of caches
    model_draws <- rbinom(n_sims, size = 1, prob = probs[2])
  } else {
    model_draws <-  sample(0:(length(all_datasets) - 1), size = n_sims, prob =  probs, replace = TRUE)
  }
  all_variables <- purrr::map(all_datasets, \(x) x$variables)
  combined_variables <- combine_draws_matrix_for_bf(all_variables,
                                                    model_draws,
                                                    NA_raw_dm = TRUE,
                                                    model_var = model_var)

  all_var_attributes <- purrr::map(all_datasets, \(x) x$var_attributes)
  combined_var_attributes <- combine_var_attributes_for_bf(all_variables,
                                                           all_var_attributes,
                                                           model_var = model_var)

  combined_generated <- list()
  for(i in 1:n_sims) {
    combined_generated[[i]] <- all_datasets[[model_draws[i] + 1]]$generated[[i]]
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
  all_submodels <- as.integer(extract_attribute_arguments_stats("submodel", results_bf$stats$attributes))
#  }
  unique_submodels <- sort(unique(all_submodels[!is.na(all_submodels)]))
  if(length(unique_submodels) <= 1) {
    stop("At least two different submodels need to be marked in the stats")
  }
  split_res <- list()
  for(submodel_id in unique_submodels) {
    split_res[[paste0("stats_H", submodel_id)]] <- get_stats_for_submodel(results_bf$stats, submodel_id)
  }
  split_res
}

# Constructs a CDF data.frame for SBC_posterior_cdf from posterior probability
# of a binary variable.
binary_to_cdf <- function(variable_name, prob1, simulated_value) {
  if(!all(simulated_value %in% c(0,1))) {
    warning("binary_to_cdf expects the simulated_value to be either 0 or 1")
  }
  discrete_to_cdf(variable_name, c(1 - prob1, prob1), simulated_value)
}

# Constructs a CDF data.frame for SBC_posterior_cdf from posterior probability
# of a discrete variable.
discrete_to_cdf <- function(variable_name, probs, simulated_value) {
  max_val <- length(probs) - 1
  if(!all(simulated_value %in% 0:max_val)) {
    warning(paste0("discrete_to_cdf expects the simulated_value to be an integer between 0 and ", max_val))
  }
  stopifnot(abs(sum(probs) - 1) < 1e-7)
  cdf_low <- c(0, cumsum(probs))
  cdf_high <- c(cumsum(probs[1:(length(probs) - 1)]), 1)
  data.frame(variable = variable_name,
             cdf_low = cdf_low[simulated_value + 1],
             cdf_high = cdf_high[simulated_value + 1])
}
