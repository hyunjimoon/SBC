#' Backend using bridgesampling to do Bayes factor calculation (and Bayesian Model Averaging) over two candidate
#' models.
#'
#' @param ... passed to [bridgesampling::bridge_sampler()].
#' @export
SBC_backend_bridgesampling <- function(..., model_var = "model", prior_probs = NULL, prior_prob1 = NULL, bridgesampling_args = list()) {

  require_package_version("bridgesampling", version = "1.0", purpose = " to use the bridgesampling SBC backend")

  all_backends <- list(...)
  if(!is.null(prior_prob1)) {
    stopifnot(length(all_backends) == 2)
    stopifnot(is.null(prior_probs))
    prior_probs <- c(1 - prior_prob1, prior_prob1)
  }

  if(is.null(prior_probs)) {
    prior_probs <- rep(1 / length(all_backends), times = length(all_backends))
  } else {
    stopifnot(is.numeric(prior_probs))
    stopifnot(length(prior_probs) == length(all_datasets))
  }

  structure(list(all_backends = all_backends,
                 model_var = model_var,
                 prior_probs = prior_probs,
                 bridgesampling_args = bridgesampling_args),
            class = "SBC_backend_bridgesampling")
}

#' Convert a fit for a given backend into a bridge sampler object.
#'
#' @description
#' Users should typically not call this method and instead let
#' [SBC_backend_bridgesampling()] to call it for them.
#'
#'
#' @param backend the underlying backend
#' @param fit corresponding to the backend
#' @param generated the generated data used to fit the underlying model
#' @param ... passed to [bridgesampling::bridge_sampler()].
#' @returns an object of class `bridge` or `bridge_list`.
#' @seealso [bridgesampling::bridge_sampler()]
#' @export
SBC_fit_to_bridge_sampler <- function(backend, fit, generated, ...) {
  UseMethod("SBC_fit_to_bridge_sampler")
}

#' @export
SBC_fit_to_bridge_sampler.default <- function(backend, fit, generated, ...) {
  all_classes <- paste0("'", class(backend), "'", collapse = ", ")
  stop(paste0("To use bridgesampling with backend of class ", all_classes, ", ",
       "you need to implement a corresponding S3 method for 'SBC_fit_to_bridge_sampler' ",
       "(most likely this would be 'SBC_fit_to_bridge_sampler.", class(backend)[1], "'"))
}

#' @export
SBC_fit_to_bridge_sampler.SBC_backend_rstan_sample <- function(backend, fit, generated, ...) {
  library(rstan)
  bridgesampling::bridge_sampler(fit, ...)
}

#' @export
SBC_fit_to_bridge_sampler.SBC_backend_brms <- function(backend, fit, generated, ...) {
  bridgesampling::bridge_sampler(fit, ...)
}

hypothesis_output_prefix <- function(i) {
  paste0("[[H",i,"]]: ")
}

hypothesis_output_prefix_bridge <- function(i) {
  paste0("[[H",i," bridge]]: ")
}


#' @export
SBC_fit.SBC_backend_bridgesampling <- function(backend, generated, cores) {
  fit_single <- function(backend_H, i) {
    # Prefix all outputs with relevant markers for post-processing
    fit_with_outputs <-
      prefix_captured(
        capture_all_outputs(SBC_fit(backend_H, generated, cores)),
        hypothesis_output_prefix(i)
      )

    reemit_captured(fit_with_outputs)
    fit <- fit_with_outputs$res


    bridge_with_outputs <-
      prefix_captured(
        capture_all_outputs(
          do.call("SBC_fit_to_bridge_sampler", c(list(backend_H, fit, generated), backend$bridgesampling_args))
        ),
        hypothesis_output_prefix_bridge(i)
      )

    list(fit = fit_with_outputs$res, bridge = bridge_with_outputs$res)
  }

  fit_bridges <- purrr::map2(backend$all_backends, 0:(length(backend$all_backends) - 1), fit_single)

  fits <- purrr::map(fit_bridges, \(x) x$fit)
  bridges <- purrr::map(fit_bridges, \(x) x$bridge)

  structure(list(
    fits = fits,
    bridges = bridges,
    model_var = backend$model_var,
    prior_probs = backend$prior_probs
  ), class = "SBC_fit_bridgesampling")
}

SBC_fit_bridgesampling_to_probs <- function(fit, log.p = FALSE) {
  # Using median-based BF when multiple bridgesampling iterations were used
  logmls <- purrr::map_dbl(fit$bridges, \(x) { median(x$logml, na.rm = TRUE)})

  if(any(is.na(logmls))) {
    print(fit$bridges)
    stop("Some logml values are NA.")
  }

  prior_log <- log(fit$prior_probs)
  log_probs_rel <- logmls + prior_log

  # Softmax to probs
  log_probs <- log_probs_rel - log_sum_exp(log_probs_rel)
  if(log.p) {
    log_probs
  } else {
    exp(log_probs)
  }
}

#' @export
SBC_fit_to_draws_matrix.SBC_fit_bridgesampling <- function(fit) {
  all_dms <- purrr::map(fit$fits, SBC_fit_to_draws_matrix)
  all_draws <- purrr::map(all_dms, posterior::merge_chains)

  all_ndraws <- purrr::map_int(all_draws, posterior::ndraws)
  shared_ndraws <- unique(all_ndraws)

  if(length(shared_ndraws) > 1) {
    warning("Unequal number of draws for each bridgesampling fit. Will subset to the smaller number.")
    shared_ndraws <- min(shared_ndraws)
    all_draws <- purrr::map(all_draws, \(x) {
      posterior::subset_draws(x, draw = 1:shared_ndraws)
    })
  }

  probs <- SBC_fit_bridgesampling_to_probs(fit)

  if(length(all_draws) == 2) {
    # Keeping the old way for 2 models to not invalidate older results
    model_draws <- rbinom(n = shared_ndraws, size = 1, prob = probs[2])
  } else {
    model_draws <- sample(0:(length(all_draws) - 1), size  = shared_ndraws, prob = probs, replace = TRUE)
  }

  combined_draws <- combine_draws_matrix_for_bf(all_draws, model_draws, model_var = fit$model_var)

  return(combined_draws)
}

#' @export
SBC_fit_specific_dquants.SBC_fit_bridgesampling <- function(fit) {
  probs <- SBC_fit_bridgesampling_to_probs(fit)
  max_index <- which.max(probs)
  is_var_name <- paste0("is_", fit$model_var, max_index - 1)
  top_var_name <- paste0("top_", fit$model_var)

  dq_args <- list(rlang::parse_quo(is_var_name, rlang::current_env()))
  names(dq_args) <- top_var_name
  dq_args$.var_attributes <- var_attributes_from_list(top_var_name, list(c(
    binary_var_attribute(), possibly_constant_var_attribute()
  )))
  do.call(derived_quantities, dq_args)

}

#' @export
SBC_posterior_cdf.SBC_fit_bridgesampling <- function(fit, variables) {
  if(fit$model_var %in% names(variables)) {
    probs <- SBC_fit_bridgesampling_to_probs(fit)
    model_cdf <- discrete_to_cdf(fit$model_var, probs, variables[fit$model_var])
    if(length(probs) == 2) {
      return(model_cdf)
    } else {
      is_model_cdf_list <- list()
      for(i in 1:length(probs)) {
        is_var_name <- paste0("is_", fit$model_var, i - 1)
        if(is_var_name %in% names(variables)) {
          is_model_cdf_list[[i]] <- binary_to_cdf(is_var_name, probs[i], variables[is_var_name])
        } else {
          is_model_cdf_list[[i]] <- NULL
        }
      }
      is_model_cdf <- do.call(rbind, is_model_cdf_list)

      max_index <- which.max(probs)
      simulated_value_top <- variables[fit$model_var] == max_index - 1
      top_prediction_cdf <- binary_to_cdf(paste0("top_", fit$model_var), probs[max_index], simulated_value_top)

      return(rbind(model_cdf, top_prediction_cdf, is_model_cdf))
    }
  } else {
    return(NULL)
  }
}

#' @export
SBC_fit_to_diagnostics.SBC_fit_bridgesampling <- function(fit, fit_output, fit_messages, fit_warnings) {

  get_prefixed_lines <- function(prefix, lines) {
    if(is.null(lines)) {
      character()
    } else {
      with_prefix <- lines[startsWith(lines, prefix)]
      prefix_removed <- substring(with_prefix, first = length(prefix) + 1)
      prefix_removed
    }
  }

  process_diag_single <- function(fit, model_index) {
    diags <- SBC_fit_to_diagnostics(fit,
                           get_prefixed_lines(hypothesis_output_prefix(model_index), fit_output),
                           get_prefixed_lines(hypothesis_output_prefix(model_index), fit_messages),
                           get_prefixed_lines(hypothesis_output_prefix(model_index), fit_warnings))
    if(!is.null(diags)) {
      names(diags) <- paste0(names(diags), "_H", model_index)
    }
    diags
  }

  model_indices <- 0:(length(fit$fits) - 1)

  diags_all <- purrr::map2(fit$fits, model_indices, process_diag_single)

  probs <- SBC_fit_bridgesampling_to_probs(fit)
  log_probs <- SBC_fit_bridgesampling_to_probs(fit, log.p = TRUE)

  get_percentage_error <- function(bridge) {
    errm <- bridgesampling::error_measures(bridge)
    if(inherits(bridge, "bridge_list")) {
      base <- bridgesampling::logml(bridge)
      perc_raw <- errm$IQR / base
      return(perc_raw * 100)
    } else {
      if(grepl("^[0-9]*%$", errm$percentage)) {
        return(as.numeric(gsub("%", "", errm$percentage)))
      } else {
        warning("Urecognized percentage format: ", errm$percentage)
        return(NA_real_)
      }
    }
  }

  get_bridge_diagnostics <- function(bridge, i) {
    bridge_warnings <- get_prefixed_lines(hypothesis_output_prefix_bridge(i), fit_warnings)

    # Example warning:
    # "209 of the 30000 log_prob() evaluations on the proposal draws produced -Inf/Inf."
    inf_proposal_lines <- grepl("evaluations.*proposal.*Inf", bridge_warnings)
    if(any(inf_proposal_lines)) {
      if(sum(inf_proposal_lines) > 1) {
        warning("Multiple lines with infinite proposals warning.")
      }
      inf_proposal_text <- (bridge_warnings[inf_proposal_lines])[1]
      n_inf_proposals <- as.integer(sub("\\D*(\\d+).*", "\\1", inf_proposal_text))
      if(n_inf_proposals <= 0) {
        warning("Problems parsing infinite proposals warnings")
      }
    } else {
      n_inf_proposals <- 0
    }
    diags <- data.frame(
      bs_error = get_percentage_error(bridge),
      n_inf_proposals = n_inf_proposals,
      logml_maxiter = any(grepl("logml.*maxiter", bridge_warnings)) #logml could not be estimated within maxiter
    )
    names(diags) <- paste0(names(diags), "_H", i)
    diags
  }

  bridge_diags_all <- purrr::map2(fit$bridges, model_indices, get_bridge_diagnostics)

  if(length(fit$fits) == 2) {
    probs_df <- data.frame(prob_H1 = probs[2], log_prob_H1 = log_probs[2])
  } else {
    probs_diag_vec <- c(probs, log_probs)
    names(probs_diag_vec) <- c(paste0("prob_H", model_indices), paste0("log_prob_H", model_indices))
    probs_df <- t(as.data.frame(probs_diag_vec))
    rownames(probs_df) <- NULL
  }
  diags_bs <- do.call(cbind, c(list(probs_df),
                               bridge_diags_all, diags_all))

  class(diags_bs) <- c("SBC_bridgesampling_diagnostics", class(diags_bs))
  attr(diags_bs, "submodel_classes") <- purrr::map(diags_all, class)
  return(diags_bs)
}

#' @export
diagnostic_types.SBC_bridgesampling_diagnostics <- function(diags) {
  submodel_classes <- attr(diags, "submodel_classes", exact = TRUE)
  if(is.null(submodel_classes) ||
     !is.list(submodel_classes)) {
    warning(
r"(The 'submodel_classes' attribute of an SBC_bridgesampling_diagnostics data.frame
is not set or is in incorrect format.
Maybe you have modified the $backend_diagnostics element of SBC_results?
If not, please file an issue at https://github.com/hyunjimoon/SBC/issues/
)")
    submodel_diags_all <- list()
    n_submodels <- 2 #Typical guess, won't hurt reporting
  } else {
    get_submodel_diags <- function(i) {
      bs_specific_diags <- list()
      bs_specific_diags[[paste0("bs_error_H", i)]] <-
        numeric_diagnostic(paste0("relative error of marginal likelihood for H", i), report = "max", error_above = 5, unit = "%")
      bs_specific_diags[[paste0("n_inf_proposals_H", i)]] <-
        count_diagnostic(paste0("H", i, ": log_prob() evaluations on the proposal draws produced -Inf/Inf"), error_above = 0, error_only = TRUE)
      bs_specific_diags[[paste0("logml_maxiter_H", i)]] <-
        logical_diagnostic(ok_value = FALSE, true_label = paste0("H", i, ": logml could not be estimated within maxiter"),
                           error_only = TRUE)

      H_diags_selected <-
          dplyr::select(diags, tidyselect::ends_with(paste0("_H",i)) & !tidyselect::all_of(c("prob_H1", "log_prob_H1", names(bs_specific_diags))))

      H_diags <- dplyr::rename_with(
        H_diags_selected,
        \(name) gsub(paste0("_H",i,"$"), "", name))


      class(H_diags) <- submodel_classes[[i + 1]]
      types_sub <- diagnostic_types(H_diags)
      types_mapped <- purrr::map(types_sub,
                                 \(diag) submodel_diagnostic(paste0("H", i), diag))
      names(types_mapped) <- paste0(names(types_sub), "_H", i)
      c(
        types_mapped,
        bs_specific_diags
      )
    }

    n_submodels <- length(submodel_classes)
    submodel_diags_all <- do.call(c,
                                  purrr::map(0:(n_submodels - 1), get_submodel_diags))
  }

  if(n_submodels == 2) {
    prob_diags <- list(
      prob_H1 = numeric_diagnostic("posterior probability of H1", report = "quantiles"),
      log_prob_H1 = skip_diagnostic()
    )
  } else {
    prob_diags <- list()
    for(i in 0:(n_submodels - 1)) {
      prob_diags[[paste0("prob_H", i)]] <- numeric_diagnostic(paste0("posterior probability of H",i), report = "quantiles")
      prob_diags[[paste0("log_prob_H", i)]] <- skip_diagnostic()
    }
  }

  c(prob_diags, submodel_diags_all)
}

#' Custom rbind implementation maintainig information about submodels
#' @exportS3Method base::rbind
rbind.SBC_bridgesampling_diagnostics <- function(...) {

  args <- list(...)

  # Working around the special dispatch for rbind
  args_class_removed <- purrr::map(args, \(x) {
    if(is.null(x)) {
      NULL
    } else {
      class(x) <- setdiff(class(x), "SBC_bridgesampling_diagnostics")
      x
    }
  })
  res <- do.call(rbind, args_class_removed)
  class(res) <- c("SBC_bridgesampling_diagnostics", class(res))

  submodel_classes_list <- purrr::map(args, \(x) attr(x, "submodel_classes", exact = TRUE))
  unique_submodel_classes <- unique(submodel_classes_list)
  if(length(unique_submodel_classes) > 1) {
    warning("Non-unique submodel classes when binding diagnostics")
  }
  submodel_classes <- unique_submodel_classes[[1]]
  if(!is.null(submodel_classes)) {
    attr(res, "submodel_classes") <- submodel_classes
  }
  res
}

#' Custom select implementation maintaining information about submodels
#' @exportS3Method dplyr::select
select.SBC_bridgesampling_diagnostics <- function(diags, ...) {
  selected <- NextMethod()
  attr(selected, "submodel_classes") <- attr(diags, "submodel_classes")
  selected
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_bridgesampling <- function(backend) {
  backend_for_hash <- backend
  backend_for_hash$all_backends <- purrr::map(backend$all_backends, SBC_backend_hash_for_cache)
  rlang::hash(backend_for_hash)
}
