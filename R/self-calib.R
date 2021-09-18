##' Auto calibrate the initial prior samples using SBC iteration
##'
##' @param generator function that generates datasets given each value in `param`
##' @param clamp_dist scalar clamp_disteter values e.g. prior_width
##' @param param rvars<S>[N] prior values i.e. parameter values to be tested "draws_rvars"
##' @param predictor rvars<N>[S] predictor values
##' @param backend stan_model that samples posterior with `datasets` simulated from the former two
##' @param target_vars function of parameters with which SBC iteration convergence are judged
##' @param scm empirical distance between the previous and the current parameter values, default is `cjs_dist`
##' @param thin scalar select 1 out of `thin` number of samples, `iter_sampling`=`thin`*M, default: 10
##' @param cnt function of parameters with which SBC iteration convergence are judged
##' @param evolve_df datafraame holding `median, mad (or sd)` of samples from every iteration
##' @param delivDir save location of output result
##' @return  next param values summarized from `S` * `M` posterior samples
##' @export

self_calib <- function(generator, clamp_val, clamp_dist, param, predictor, backend, target_vars, thin, cnt, evolve_df, delivDir, type = type){
  # adjust draws (='S') of clamp_dist and predictor to that of new param
  S <- niterations(param[[1]])
  clamp_dist <- subset_draws(clamp_dist, draw = 1:S)
  predictor <- subset_draws(predictor, draw = 1:S)
  message(paste("self_calib iter", cnt))
  message(paste("S", S))
  # generate-inference p_post(theta) = f(theta'|y) * p(y|theta)
  result <- compute_results(generator(clamp_val, clamp_dist, param, predictor), backend, thin_ranks = thin)
  # update
  param_next <- update_param(param, result, thin, cnt, delivDir, type = type)
  # save
  summ <- summarise_draws(param, median, sd) %>% filter(variable == target_vars)
  for (tv in target_vars){
    evolve_df[[tv]]$median[cnt] <- filter(summ, variable == tv)["median"]
    evolve_df[[tv]]$sd[cnt] <-  filter(summ, variable == tv)["sd"]
    evolve_df[[tv]]$scm <- cjs_dist(param[[tv]], param_next[[tv]])
    print(paste(tv, evolve_df[[tv]]$median[cnt]))
    if(cnt == 0){evolve_df[[tv]]$scm_init <- evolve_df[[tv]]$scm}
  }
  # terminate
  if(iter_stop(param, param_next, target_vars, lapply(evolve_df, '[[', 'scm'))){ # S-S > S-4000S (stable compare)
    csv_save(evolve_df, delivDir, type = "evolve")
    return (param_next)
  }
  cnt = cnt + 1
  print(param_next)
  return(self_calib(generator, clamp_val, clamp_dist, param_next, predictor, backend, target_vars, thin, cnt, evolve_df, delivDir, type = type))
}

##' Judge whether the SBC iteration have converged
##'
##' @param param numeric vector of parameter values to be tested a.k.a prior sample
##' @param param_next numeric vector of updated parameter values after one SBC cycle a.k.a posterior values
##' @param target_vars function of parameters with which SBC iteration convergence are judged, null is strict testing all parameters
##' @param scm_init initial distance between initial parameter and its proposal
##' @return boolean, judge convergence based on scm-scm_init ratio (relative) or `param`-`param_next` distance (absolute)
##' @export
iter_stop <- function(param, param_next, target_vars, scm_init){
      return(all(unlist(lapply(target_vars,
                               FUN = function(tv) cjs_dist(param[[tv]], param_next[[tv]]) < 0.5 * unlist(scm_init[[tv]])))))
}

# Updating parameter value `param` to `param_next` for valid comparison and to control `S`
# Resample `param_next` sample matrix for with PIT weight of F_{prior}F_{post}^{-1}
# S * n_sample posterior S as comparison threshold is possible for the same number of samples
##'
##' @param param numeric vector of prior values i.e. parameter values to be tested
##' @param post posterior samples of length S * M

##' @return resampled posterior with prior information
##' @export
update_param <-function(theta_, result, thin, cnt, delivDir, target_vars = names(theta_), type = "kde", kde_bandwidth=0.5){
  S <- niterations(theta_[[1]])
  tf <- list()
  post <- list()
  theta <- list()
  theta_tf <- list()
  theta_tf_next <- list()
  for (tv in target_vars){
    agg_post <- function(s){
      if(!is.null(result$fits[[s]])){# ADVI fit can be NULL
        return (thin_draws(subset_draws(SBC_fit_to_draws_matrix(result$fits[[s]]), variable = tv), thin))}
    }
    for (s in 1:S){
      if(tv == target_vars[1] & !is.null(result$fits[[s]])){
        M <- length(agg_post(s))
      }
    }
    theta[[tv]] <- unlist(lapply(seq(1: S), agg_post)) # can reject/thin unideal nrow = S, ncol = M
  }
  S <- floor(length(theta[[1]])/M)
  pp_overlay_save(theta_, as_draws_rvars(theta), cnt, delivDir)
  out <- tf_param(theta)
  theta_tf <-  tibble::as_tibble(out$param)
  tf <- out$tf
  # parameteric
  if(type == "kde"){
    sample_kde <- function(x, bandwidth, n_samples){
      kde <- density(x, bw=bandwidth)
      bw <- kde$bw
      draws <- sample(x, n_samples, replace = TRUE)
      rnorm(n_samples, draws, bw) # sample(x) + normal(0, bandwidth)
    }
    for(tv in target_vars){
      theta_tf_next[[tv]] <- sample_kde(dplyr::pull(theta_tf, tv), bandwidth = kde_bandwidth, n_samples = S)
    }
  }else if (type == "sample"){
    for(tv in target_vars){
      theta_tf_next[[tv]] <- theta_tf[[tv]][unlist(lapply(c(1:S), FUN = function(i){sample(size = 1, 1:M) + (i - 1) * M}))]
    }
    # nonparameteric
  }else if (type == "bin"){
    # to be implemented
    theta_tf <- theta_tf
  }else if(type == "filter"){
    # F_{post}^{-1}*F_{prior}(param) #TODO loo_subsample.R
    for(tv in target_vars){
      datagrid <- sort(draws_of(theta_[[tv]]))
      theta_tf_next[[tv]] <- draws_of(resample_draws(as_draws(rvar(datagrid)),
                                          tabulate(ecdf(datagrid)(theta[[tv]] ) * S, nbins = S))[[1]])
    }
  }
  theta <- invtf_param(theta_tf_next, tf)
  return (as_draws_rvars(theta))
}

ar_param <- function(target_ftn = cjs_dist, prop, param, evolve_df, cnt){
  param_next <- param
  scm_prop <- target_ftn(prop, param[[tv]])
  if(cnt == 0){
    evolve_df[[tv]]$scm_init <- scm_prop
    evolve_df[[tv]]$scm <- scm_prop
  }
  scm_t <- evolve_df[[tv]]$scm
  # reject-acceptance
  if(runif(1) < log(scm_prop)/log(scm_t)){
    param_next[[tv]] <- prop
    evolve_df[[tv]]$scm <- scm_prop
  }
  return(param_next)
}

