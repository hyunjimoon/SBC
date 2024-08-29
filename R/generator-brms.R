#' Create a brms generator.
#'
#' Brms generator uses a brms model with `sample_prior = "only"` to generate
#' new datasets.
#'
#' @param formula,data,... arguments passed to `brms::brm`
#' @param generate_lp whether to compute the overall log-likelihood of the model
#' as an additional variable. This can be somewhat computationally expensive,
#' but improves sensitivity of the SBC process.
#' @param generate_lp_chunksize Will determine whether
#' the chunk size used with [future.apply::future_mapply()]. Set quite high
#' by default as the parallelism only benefits for large individual datasets/number of
#' simulations.
#' @param out_stan_file A filename for the generated Stan code. Useful for
#'    debugging and for avoiding unnecessary recompilations.
#' @export
SBC_generator_brms <- function(formula, data, ...,  generate_lp = TRUE,
                               generate_lp_chunksize = 5000 / nrow(data),
                               out_stan_file = NULL) {
  require_brms_version("brms generator")

  model_data <- brms::make_standata(formula = formula, data = data, ..., sample_prior = "only")
  class(model_data) <- NULL

  args <- list(...)

  if(!is.null(args$algorithm) && args$algorithm != "sampling" && args$algorithm != "meanfield") {
    stop("Algorithms other than sampling and meanfield not supported yet")
  }

  compiled_model <- stanmodel_for_brms(formula = formula, data = data, out_stan_file = out_stan_file, ...)



  structure(list(
    model_data = model_data,
    args = c(list(formula = formula, data = data), args),
    generate_lp = generate_lp,
    generate_lp_chunksize = generate_lp_chunksize,
    compiled_model = compiled_model),
    class = "SBC_generator_brms")
}

#' @export
generate_datasets.SBC_generator_brms <- function(generator, n_sims, n_datasets = NULL) {
  if(!is.null(n_datasets)) {
    warning("n_datasets argument is deprecated, use n_sims instead")
    if(missing(n_sims)) {
      n_sims <- n_datasets
    }
  }

  #TODO pass args for control, warmup, .... to sampling
  if(inherits(generator$compiled_model, "CmdStanModel")) {
      args_for_fitting <- translate_rstan_args_to_cmdstan(generator$args, include_unrecognized = FALSE)
      args_for_fitting$data <- generator$model_data

      if(is.null(args_for_fitting$chains)) {
        args_for_fitting$chains <- 1
      }
      if(is.null(args_for_fitting$thin)) {
        args_for_fitting$thin <- 1
      }


      args_for_fitting$iter_sampling <- ceiling(n_sims / args_for_fitting$chains) * args_for_fitting$thin

      args_for_fitting
      prior_fit <- do.call(generator$compiled_model$sample,
                           args_for_fitting)


      # Change once https://github.com/stan-dev/cmdstanr/issues/205
      # is resolved
      summ <- prior_fit$summary() # Can trigger warnings for treedepth/divergent ...
      if(any(is.na(summ$rhat)) || any(is.na(summ$ess_bulk))) {
        message("Some rhats / bulk effective sample sizes are NA.\n",
                "This is not a problem if you have constant elements in your model or you generated very few datasets and/or used little thinning.")
      }
      max_rhat <- max(c(-Inf, summ$rhat), na.rm = TRUE)
      if(max_rhat > 1.01) {
        message("Warning: Some rhats are > 1.01 indicating the prior was not explored well.\n",
                "The highest rhat is ", round(max_rhat, 2)," for ", summ$variable[which.max(summ$rhat)],
                "\nConsider adding warmup iterations (via 'warmup' argument).")
      }
      min_ess <- min(c(Inf, summ$ess_bulk), na.rm = TRUE)
      if(min_ess < n_sims / 2) {
        message("Warning: Bulk effective sample size for some parameters is less than half the number of simulations.\n",
                "The lowest ESS_bulk/n_sims is ", round(min_ess / n_sims, 2)," for ", summ$parameter[which.min(summ$ess_bulk)],
                "\nConsider increased thinning  (via 'thin' argument) .")
      }

  } else if (inherits(generator$compiled_model, "stanmodel")) {
    args_to_pass <- c("thin", "warmup", "control", "refresh")
    args_for_fitting <- c(
      list(object = generator$compiled_model,
           data = generator$model_data,
           chains = 1
      ),
      generator$args[intersect(args_to_pass, names(generator$args))]
    )

    if(is.null(args_for_fitting$warmup)) {
      args_for_fitting$warmup <- 1000
    }
    if(is.null(args_for_fitting$chains)) {
      args_for_fitting$chains <- 1
    }
    if(is.null(args_for_fitting$thin)) {
      args_for_fitting$thin <- 1
    }

    args_for_fitting$iter <- args_for_fitting$warmup + ceiling(n_sims / args_for_fitting$chains) * args_for_fitting$thin

    prior_fit <- do.call(rstan::sampling, args_for_fitting)

  } else {
    stop("Invalid generator$compiled_model")
  }

  prior_fit_brms <- brmsfit_from_stanfit(prior_fit, generator$args)

  generated <- brms_full_ppred(prior_fit_brms)

  if(generator$generate_lp) {
    ## Compute the likelihoods (observation model)
    if(generator$generate_lp_chunksize < n_sims) {
      log_likelihoods <- future.apply::future_mapply(
        FUN = function(new_dataset, i) {
          ll <- brms::log_lik(prior_fit_brms, newdata = new_dataset, draw_ids = i, cores = 1)
          sum(ll)
          },
        generated, 1:n_sims,
        future.chunk.size = generator$generate_lp_chunksize)
    } else {
      log_likelihoods <- numeric(n_sims)
      for(i in 1:n_sims) {
        ll <- log_lik(prior_fit_brms, newdata = generated[[i]], draw_ids = i, cores = 1)
        log_likelihoods[i] <- sum(ll)
      }
    }
  }

  # Add the observation likelihoods to the lp__ of the prior to
  # (hopefully) get total likelihood
  draws <- posterior::as_draws_matrix(prior_fit_brms$fit)
  if(generator$generate_lp) {
    draws[,"lp__"] <- draws[,"lp__"] + log_likelihoods
  } else {
    new_variables <- setdiff(posterior::variables(draws), "lp__")
    draws <- posterior::subset_draws(draws, variable = new_variables)
  }

  SBC_datasets(draws, generated)

}


#' Full forward sampling of a the response of brms fit, including multivariate models.
#'
#' @param fit An object of class `brmsfit`
#' @param newdata An optional data.frame for which to evaluate predictions. If NULL (default), the original data of the model is used.
#' @param draws An integer vector specifying the posterior draws to be used. If NULL (the default), all draws are used.
#' @param validate_all if TRUE, validation of input data will be done in all iterations, otherwise only once
#'
#' @return A list of data.frames containing the draws.
#'
#' @keywords internal
#' @examples # Pending
brms_full_ppred <- function(fit, newdata = NULL, draws = NULL, validate_all = FALSE) {
  # 1. determine term hierarchy
  resp <- brms_response_sequence(fit)
  # 2.1. initialize dataframe using the original fit's data
  if(is.null(newdata)) newdata <- fit$data
  n <- nrow(newdata)
  # 2.3. if no draws set, range from 1 to all iters (check draws < iters)
	if(is.null(draws)) draws <- seq_len(posterior::ndraws(fit))
  # 2.4. create list to hold data
  pp_data <- list()

  if(!validate_all) {
    # Validate once
    newdata <- brms::validate_newdata(newdata, fit)
  }

  for (i in draws) {
    pp_data[[i]] <- newdata
    for (vars in resp) {
      pp_data[[i]][, vars] <- array(
        brms::posterior_predict(
          fit, newdata = pp_data[[i]],
          resp = vars, draw_ids = i,
          skip_validate = !validate_all),
        dim = c(1, n, length(vars)))[1,,]
    }
  }
  pp_data
}

nodes_by_depth <- function(adj_matrix) {
  depth_list <- list()
  var_names <- rownames(adj_matrix)
  while(nrow(adj_matrix)) {
    pos <- which(apply(adj_matrix, 1, sum) == 0)
    depth_list <- c(depth_list, list(var_names[pos]))

    var_names <- var_names[-pos]
    adj_matrix <- adj_matrix[-pos, -pos, drop = FALSE]
  }
  depth_list
}

#' Determine the response sequence of brms model
#' @keywords internal
brms_response_sequence <- function(x){ UseMethod("brms_response_sequence") }

#' @method brms_response_sequence brmsfit
#' @export
brms_response_sequence.brmsfit <- function(x){ brms_response_sequence(x$formula) }

#' @method brms_response_sequence bform
#' @export
brms_response_sequence.bform <- function(x){
  term_list <- brms_response_sequence(brms::brmsterms(x))
  resp_vars <- names(term_list)

  adjacency <- t(sapply(term_list, \(x)is.element(resp_vars, x)))
  attr(adjacency, "dimnames") <- list(resp_vars, resp_vars)
  nodes_by_depth(adjacency)
}

#' @method brms_response_sequence mvbrmsterms
#' @export
brms_response_sequence.mvbrmsterms <- function(x){
  names(x$terms) <- NULL
  sapply(x$terms, brms_response_sequence)
}

#' @method brms_response_sequence brmsterms
#' @export
brms_response_sequence.brmsterms <- function(x){
  vars <- list(unique(unlist(lapply(x$dpars, brms_response_sequence))))
  names(vars) <- all.vars(x$respform)
  vars
}

#' @method brms_response_sequence btl
#' @export
brms_response_sequence.btl <- function(x){ c("1", all.vars(x$formula)) }
