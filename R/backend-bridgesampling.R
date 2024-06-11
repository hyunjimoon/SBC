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
  for(v in all_variables_names) {
    if(v %in% posterior::variables(dm0)) {
      var0 <- as.numeric(dm0[, v])
    } else {
      var0 <- NA_real_
    }
    if(v %in% posterior::variables(dm1)) {
      var1 <- as.numeric(dm1[, v])
    } else {
      var1 <- NA_real_
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
combine_var_attributes_for_bf <- function(dm0, dm1, var_attributes0, var_attributes1, model_var = "model") {
  stopifnot(posterior::is_draws_matrix(dm0))
  stopifnot(posterior::is_draws_matrix(dm1))

  var_attributes0 <- validate_var_attributes(var_attributes0)
  var_attributes1 <- validate_var_attributes(var_attributes1)

  stopifnot(is.character(model_var) & length(model_var) == 1)


  warning("Temporary implementation of combine_var_attributes_for_bf")
  raw_attr_vec <- c(hidden_var_attribute(), na_valid_var_attribute())

  attr0_names <- unique(variable_names_to_var_attributes_names(posterior::variables(dm0)))
  attr0_raw <- rep(list(raw_attr_vec), length(attr0_names))
  names(attr0_raw) <- paste0(".m0.",attr0_names)


  attr1_names <- unique(variable_names_to_var_attributes_names(posterior::variables(dm1)))
  attr1_raw <- rep(list(raw_attr_vec), length(attr1_names))
  names(attr1_raw) <- paste0(".m1.",attr1_names)


  attr_model <- list(c(binary_var_attribute(), possibly_constant_var_attribute()))
  names(attr_model) <- model_var

  return(c(var_attributes0, var_attributes1, attr_model, attr0_raw, attr1_raw))
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

#' Backend using bridgesampling to do Bayes factor calculation (and Bayesian Model Averaging) over two candidate
#' models.
#'
#' @param ... passed to [bridgesampling::bridge_sampler()].
#' @export
SBC_backend_bridgesampling <- function(backend_H0, backend_H1, model_var = "model", ...) {
  require_package_version("bridgesampling", version = "1.0", purpose = " to use the bridgesampling SBC backend")
  structure(list(backend_H0 = backend_H0,
                 backend_H1 = backend_H1,
                 model_var = model_var,
                 bridgesampling_args = list(...)),
            class = "SBC_backend_bridgesampling")
}

#' Convert a fit for a given backend into a bridge sampler object.
#' @param ... passed to [bridgesampling::bridge_sampler()].
#' @seealso [bridgesampling::bridge_sampler()]
#' @export
SBC_fit_to_bridge_sampler <- function(backend, fit, generated, ...) {
  UseMethod("SBC_fit_to_bridge_sampler")
}

#' @export
SBC_fit_to_bridge_sampler.default <- function(backend, fit, generated, ...) {
  stop(paste0("To use bridgesampling with backend of class '", class(backend)), "',",
       "you need to implement a corresponding S3 method for SBC_fit_to_bridge_sampler")
}

#' @export
SBC_fit_to_bridge_sampler.SBC_backend_rstan_sample <- function(backend, fit, generated, ...) {
  bridgesampling::bridge_sampler(fit, ...)
}

#' @export
SBC_fit_to_bridge_sampler.SBC_backend_brms <- function(backend, fit, ...) {
  bridgesampling::bridge_sampler(fit, ...)
}


#' @export
SBC_fit.SBC_backend_bridgesampling <- function(backend, generated, cores) {
  fit0 <- SBC_fit(backend$backend_H0, generated, cores)
  fit1 <- SBC_fit(backend$backend_H1, generated, cores)
  bridge_H0 <- do.call("SBC_fit_to_bridge_sampler", c(list(backend$backend_H0, fit0, generated), backend$bridgesampling_args))
  bridge_H1 <- do.call("SBC_fit_to_bridge_sampler", c(list(backend$backend_H1, fit1, generated), backend$bridgesampling_args))
  structure(list(
    fit0 = fit0,
    fit1 = fit1,
    bridge_H0 = bridge_H0,
    bridge_H1 = bridge_H1
  ), class = "SBC_fit_bridgesampling")
}

#' @export
SBC_fit_to_draws_matrix.SBC_fit_bridgesampling <- function(fit) {
  draws0 <- posterior::merge_chains(SBC_fit_to_draws_matrix(fit$fit0))
  draws1 <- posterior::merge_chains(SBC_fit_to_draws_matrix(fit$fit1))

  if(posterior::ndraws(draws0) != posterior::ndraws(draws1)) {
    warning("Unequal number of draws for each bridgesampling fit. Will subset to the smaller number.")
    if(posterior::ndraws(draws0) > posterior::ndraws(draws1)) {
      draws0 <- posterior::subset_draws(draws0, draw = 1:posterior::ndraws(draws1))
    } else {
      draws1 <- posterior::subset_draws(draws1, draw = 1:posterior::ndraws(draws0))
    }
  }

  logml0 <- bridgesampling::logml(fit$bridge_H0)
  logml1 <- bridgesampling::logml(fit$bridge_H1)


  log_bf_01 <- bridgesampling::bf(fit$bridge_H0, fit$bridge_H1, log = TRUE)$bf
  prob1 <- plogis(-log_bf_01)

  total_draws <- posterior::ndraws(draws0)

  model_draws <- rbinom(n = total_draws, size = 1, prob = prob1)

  combined_draws <- combine_draws_matrix_for_bf(draws0, draws1, model_draws)

  return(combined_draws)
}

#' @export
SBC_fit_to_diagnostics.SBC_fit_bridgesampling <- function(fit, fit_output, fit_messages, fit_warnings) {
  diags0 <- SBC_fit_to_diagnostics(fit$fit0, fit_output, fit_messages, fit_warnings)
  diags1 <- SBC_fit_to_diagnostics(fit$fit1, fit_output, fit_messages, fit_warnings)


  log_bf_01 <- bridgesampling::bf(fit$bridge_H0, fit$bridge_H1, log = TRUE)$bf
  prob1 <- plogis(-log_bf_01)

  percentage_error0 <- bridgesampling::error_measures(fit$bridge_H0)$percentage
  percentage_error1 <- bridgesampling::error_measures(fit$bridge_H1)$percentage
  diags_bs <- data.frame(prob_H1 = prob1, bs_error_H0 = percentage_error0, bs_error_H1 = percentage_error1)

  if(!is.null(diags0)) {
    names(diags0) <- paste0(names(diags0), "_H0")
    diags_bs <- cbind(diags_bs, diags0)
  }

  if(!is.null(diags1)) {
    names(diags1) <- paste0(names(diags1), "_H1")
    diags_bs <- cbind(diags_bs, diags1)
  }

  return(diags_bs)
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_bridgesampling <- function(backend) {
  backend_for_hash <- backend
  backend_for_hash$backend_H0 <- SBC_backend_hash_for_cache(backend$backend_H0)
  backend_for_hash$backend_H1 <- SBC_backend_hash_for_cache(backend$backend_H1)
  rlang::hash(backend_for_hash)
}
