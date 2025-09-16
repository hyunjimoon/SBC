#' Backend using bridgesampling to do Bayes factor calculation (and Bayesian Model Averaging) over two candidate
#' models.
#'
#' @param ... passed to [bridgesampling::bridge_sampler()].
#' @export
SBC_backend_bridgesampling <- function(backend_H0, backend_H1, model_var = "model", prior_prob1 = 0.5, ...) {
  require_package_version("bridgesampling", version = "1.0", purpose = " to use the bridgesampling SBC backend")
  structure(list(backend_H0 = backend_H0,
                 backend_H1 = backend_H1,
                 model_var = model_var,
                 prior_prob1 = prior_prob1,
                 bridgesampling_args = list(...)),
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

  fit_bridge_0 <- fit_single(backend$backend_H0, 0)
  fit_bridge_1 <- fit_single(backend$backend_H1, 1)

  structure(list(
    fit0 = fit_bridge_0$fit,
    fit1 = fit_bridge_1$fit,
    bridge_H0 = fit_bridge_0$bridge,
    bridge_H1 = fit_bridge_1$bridge,
    model_var = backend$model_var,
    prior_prob1 = backend$prior_prob1
  ), class = "SBC_fit_bridgesampling")
}

SBC_fit_bridgesampling_to_prob1 <- function(fit, log.p = FALSE) {
  bf_res <- bridgesampling::bf(fit$bridge_H0, fit$bridge_H1, log = TRUE)
  if(inherits(bf_res, "bf_bridge_list")) {
    log_bf_01 <- bf_res$bf_median_based
  } else {
    log_bf_01 <- bf_res$bf
  }
  if(is.na(log_bf_01)) {
    print(fit$bridge_H0)
    print(fit$bridge_H1)
    stop("Bayes factor is NA.")
  }
  prior_log <- log(fit$prior_prob1) - log1p( -fit$prior_prob1)
  prob1 <- plogis(-log_bf_01 + prior_log, log.p = log.p)
  return(prob1)
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

  prob1 <- SBC_fit_bridgesampling_to_prob1(fit)

  total_draws <- posterior::ndraws(draws0)

  model_draws <- rbinom(n = total_draws, size = 1, prob = prob1)

  combined_draws <- combine_draws_matrix_for_bf(draws0, draws1, model_draws, model_var = fit$model_var)

  return(combined_draws)
}

#' @export
SBC_posterior_cdf.SBC_fit_bridgesampling <- function(fit, variables) {
  if(fit$model_var %in% names(variables)) {
    prob1 <- SBC_fit_bridgesampling_to_prob1(fit)
    return(binary_to_cdf(fit$model_var, prob1, variables[fit$model_var]))
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

  diags0 <- SBC_fit_to_diagnostics(fit$fit0,
                                   get_prefixed_lines(hypothesis_output_prefix(0), fit_output),
                                   get_prefixed_lines(hypothesis_output_prefix(0), fit_messages),
                                   get_prefixed_lines(hypothesis_output_prefix(0), fit_warnings))
  diags1 <- SBC_fit_to_diagnostics(fit$fit1,
                                   get_prefixed_lines(hypothesis_output_prefix(1), fit_output),
                                   get_prefixed_lines(hypothesis_output_prefix(1), fit_messages),
                                   get_prefixed_lines(hypothesis_output_prefix(1), fit_warnings))

  prob1 <- SBC_fit_bridgesampling_to_prob1(fit)
  log_prob1 <- SBC_fit_bridgesampling_to_prob1(fit, log.p = TRUE)

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

  diags_bs <- cbind(
    data.frame(prob_H1 = prob1, log_prob_H1 = log_prob1),
    get_bridge_diagnostics(fit$bridge_H0, 0),
    get_bridge_diagnostics(fit$bridge_H1, 1)
  )

  if(!is.null(diags0)) {
    names(diags0) <- paste0(names(diags0), "_H0")
    diags_bs <- cbind(diags_bs, diags0)
  }

  if(!is.null(diags1)) {
    names(diags1) <- paste0(names(diags1), "_H1")
    diags_bs <- cbind(diags_bs, diags1)
  }

  class(diags_bs) <- c("SBC_bridgesampling_diagnostics", class(diags_bs))
  attr(diags_bs, "submodel_classes") <- list(H0 = class(diags0), H1 = class(diags1))
  return(diags_bs)
}

#' @export
diagnostic_types.SBC_bridgesampling_diagnostics <- function(diags) {
  submodel_classes <- attr(diags, "submodel_classes", exact = TRUE)
  if(is.null(submodel_classes) ||
     !is.list(submodel_classes) || !identical(names(submodel_classes), c("H0", "H1"))) {
    warning(
r"(The 'submodel_classes' attribute of an SBC_bridgesampling_diagnostics data.frame
is not set or is in incorrect format.
Maybe you have modified the $backend_diagnostics element of SBC_results?
If not, please file an issue at https://github.com/hyunjimoon/SBC/issues/
)")
    submodel_diags <- list()
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


      class(H_diags) <- submodel_classes[[paste0("H",i)]]
      types_sub <- diagnostic_types(H_diags)
      types_mapped <- purrr::map(types_sub,
                                 \(diag) submodel_diagnostic(paste0("H", i), diag))
      names(types_mapped) <- paste0(names(types_sub), "_H", i)
      c(
        types_mapped,
        bs_specific_diags
      )
    }

    submodel_diags <- c(
      get_submodel_diags(0),
      get_submodel_diags(1)
    )

  }

  c(
    list(
      prob_H1 = numeric_diagnostic("posterior probability of H1", report = "quantiles"),
      log_prob_H1 = skip_diagnostic()
    ),
    submodel_diags
  )
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

#' Custom select implementation maintainig information about submodels
#' @exportS3Method dplyr::select
select.SBC_bridgesampling_diagnostics <- function(diags, ...) {
  selected <- NextMethod()
  attr(selected, "submodel_classes") <- attr(diags, "submodel_classes")
  selected
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_bridgesampling <- function(backend) {
  backend_for_hash <- backend
  backend_for_hash$backend_H0 <- SBC_backend_hash_for_cache(backend$backend_H0)
  backend_for_hash$backend_H1 <- SBC_backend_hash_for_cache(backend$backend_H1)
  # Keep caches from older versions valid
  if(!is.null(backend$prior_prob1) && backend$prior_prob1 == 0.5) {
    backend_for_hash$prior_prob1 <- NULL
  }
  rlang::hash(backend_for_hash)
}
