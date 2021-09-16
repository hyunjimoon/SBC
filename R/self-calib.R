##' Auto calibrate the initial prior samples using SBC iteration
##'
##' @param generator function that generates datasets given each value in `param`
##' @param hyperparam scalar hyperparameter values e.g. prior_width
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

self_calib <- function(generator, hyperparam, param, predictor, backend, target_vars, thin, cnt, evolve_df, delivDir, type = "all"){
  message(paste("self_calib iter", cnt))
  S = niterations(param[[1]])
  # generate-inference p_post(theta) = f(theta'|y) * p(y|theta)
  result <- compute_results(generator(hyperparam, param, predictor), backend, thin_ranks = thin)
  # proposal
  prop <- prop_param(param, result, thin, type = type, cnt, delivDir)
  # accept-reject
  param_next <- prop #ar_param(param, prop)
  summ <- summarise_draws(param, median, sd) %>% filter(variable == target_vars)
  # save
  for (tv in target_vars){
    evolve_df[[tv]]$median[cnt] <- filter(summ, variable == tv)["median"]
    evolve_df[[tv]]$sd[cnt] <-  filter(summ, variable == tv)["sd"]
    evolve_df[[tv]]$scm <- cjs_dist(param[[tv]], param_next[[tv]])
    print(paste(tv, evolve_df[[tv]]$median[cnt]))
    if(cnt == 0){evolve_df[[tv]]$scm_init <- evolve_df[[tv]]$scm}
  }
  pp_overlay_save(param, param_next, cnt, delivDir)
  # terminate
  if(iter_stop(param, param_next, target_vars, lapply(evolve_df, '[[', 'scm'))){ # S-S > S-4000S (stable compare)
    csv_save(evolve_df, delivDir, type = "evolve")
    return (param_next)
  }
  cnt = cnt + 1
  return(self_calib(generator, hyperparam, param_next, predictor, backend, target_vars, thin, cnt, evolve_df, delivDir, type = type))
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
prop_param <-function(param, result, thin, cnt, delivDir, target_vars = names(param), type = "bin", kde_bandwidth=0.5){
  print(type)
  S <- niterations(param[[1]])
  tf <- list()
  post <- list()
  prop_tf <- list()
  for (tv in target_vars){
    agg_post <- function(s){
      if(!is.null(result$fits[[s]])){# ADVI fit can be NULL
        return (thin_draws(subset_draws(SBC_fit_to_draws_matrix(result$fits[[s]]), variable = tv), thin))}
    }
    post[[tv]] <- unlist(lapply(seq(1: S), agg_post)) # can reject/thin unideal nrow = S, ncol = M
    if(all(post[[tv]] >0)){if(all(post[[tv]] < 1)){tf[[tv]] = "logit"} else{tf[[tv]] = "log"}}
  }
  SM <- length(post[[1]])
  post <- tibble::as_tibble(post)
  if (type == "all"){
      return(as_draws_rvars(post[sample(1:SM, S),]))
  }else if (type == "bin"){
    for(tv in target_vars){
      if(is.null(tf[[tv]])){break} #2.6e+64 for median 0.9 -> cutoff .01
      if(tf[[tv]] == "log"){post_tf[[tv]] = log(post_tf[[tv]])}
      if(tf[[tv]] == "logit"){post_tf[[tv]] = lapply(post_tf[[tv]], FUN = function(p){log(p/(1-p))})}
    }

    prop <- prop_tf # template
    for(tv in target_vars){
      if(is.null(tf[[tv]])){prop[[tv]] = prop_tf[[tv]]
      }else if(tf[[tv]] == "log"){prop[[tv]] = exp(prop_tf[[tv]])
      }else if(tf[[tv]] == "logit"){prop[[tv]] = lapply(prop_tf[[tv]], FUN = function(x){exp(x)/(1+exp(x))})
      }
    }
    return(as_draws_rvars(prop))
  }else if(type == "filter"){
    # F_{post}^{-1}*F_{prior}(param) #TODO loo_subsample.R
    S <-niterations(param[[1]])
    param_ord <- sort(draws_of(param))
    return (rvar(draws_of(resample_draws(as_draws(rvar(param_ord)),
                                         tabulate(ecdf(param_ord)(post) * S, nbins = S))[[1]])))
  }
  else if(type == "kde"){
    rvar_list <- list()
    for(tv in target_vars){
      #prop_param(param, result, thin, type = "all", cnt, delivDir)
      rvar_list[[tv]] <- rvar(sample_kde(dplyr::pull(post, tv), bandwidth = kde_bandwidth, n_samples = S))
      print(as_draws_rvars(rvar_list))
      print(as_draws_rvars(rvar_list)[[tv]])
      #return(rvar(sample_kde(post)))
    }
    return(as_draws_rvars(rvar_list))
  }

}


sample_kde <- function(x, bandwidth, n_samples){
  kde <- density(x, bw=bandwidth)
  bw <- kde$bw
  draws <- sample(x, n_samples, replace = TRUE)
  rnorm(n_samples, draws, bw) # sample(x) + normal(0, bandwidth)
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

