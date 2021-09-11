##' Auto calibrate the initial prior samples using SBC iteration
##'
##' @param generator function that generates datasets given each value in `param`
##' @param hyperparam scalar hyperparameter values e.g. prior_width
##' @param param rvars<S>[N] prior values i.e. parameter values to be tested "draws_rvars"
##' @param predictor rvars<N>[S] predictor values
##' @param backend stan_model that samples posterior with `datasets` simulated from the former two
##' @param target_vars function of parameters with which SBC iteration convergence are judged
##' @param dpp empirical distance between the previous and the current parameter values, default is `cjs_dist`
##' @param thin scalar select 1 out of `thin` number of samples, `iter_sampling`=`thin`*M, default: 10
##' @param cnt function of parameters with which SBC iteration convergence are judged
##' @param evolve_df datafraame holding `median, mad (or sd)` of samples from every iteration
##' @param delivDir save location of output result
##' @return  next param values summarized from `S` * `M` posterior samples
##' @export

self_calib <- function(generator, hyperparam, param, predictor, backend, target_vars,thin, cnt, evolve_df, delivDir){
  S = niterations(param)
  # generate-inference p_post(theta) = f(theta'|y) * p(y|theta)
  result <- compute_results(generator(hyperparam, param, predictor), backend, thin = thin)
  # proposal = F_{post}^{-1}*F_{prior}(param) #TODO loo_subsample.R
  param_next <- param
  for (tv in target_vars){
    agg_post <- function(s){
      if(!is.null(result$fits[[s]])){# ADVI fit can be NULL
        return (thin_draws(subset_draws(SBC_fit_to_draws_matrix(result$fits[[s]]), variable = tv), thin))}
    }
    post_vec <- unlist(lapply(seq(1: S), agg_post)) # can reject/thin unideal nrow = S, ncol = M
    param_prop <- prop_param(param[[tv]], post_vec)
    dpp_prop <- cjs_dist(param_prop, param[[tv]])
    if(cnt == 0){
      evolve_df[[tv]]$dpp_init <- dpp_prop
      evolve_df[[tv]]$dpp <- dpp_prop
    }
    dpp_t <- evolve_df[[tv]]$dpp
    # reject-acceptance
    if(runif(1) < log(dpp_prop)/log(dpp_t)){
      param_next[[tv]] <- param_prop
      evolve_df[[tv]]$dpp <- dpp_prop
    }
    summ <- summarise_draws(param, median, sd) %>% filter(variable == tv)
    evolve_df[[tv]]$median[cnt] <- as.numeric(summ["median"])
    evolve_df[[tv]]$sd[cnt] <- as.numeric(summ["sd"])
  }
  pp_overlay_save(param, param_next, cnt)
  # terminate
  if(iter_stop(param, param_next, target_vars,  lapply(evolve_df, '[', 'dpp'))){ # S-S > S-4000S (stable compare)
    csv_save(evolve_df, delivDir, type = "evolve")
    return (param_next)
  }
  cnt = cnt + 1
  return(self_calib(generator, hyperparam, param_next, predictor, backend, target_vars,
                            thin, cnt, evolve_df, delivDir))
}

##' Judge whether the SBC iteration have converged
##'
##' @param param numeric vector of parameter values to be tested a.k.a prior sample
##' @param param_next numeric vector of updated parameter values after one SBC cycle a.k.a posterior values
##' @param target_vars function of parameters with which SBC iteration convergence are judged, null is strict testing all parameters
##' @param dpp_init initial distance between initial parameter and its proposal
##' @return boolean, judge convergence based on dpp-dpp_init ratio (relative) or `param`-`param_next` distance (absolute)
##' @export
iter_stop <- function(param, param_next, target_vars, dpp_init, type = "rel"){
  if (type == "abs"){
    post_r_loc <- lapply(param_next, mean)
    post_r_scale <- lapply(param_next, sd)
    r_loc <- list()
    r_scale <- list()
    for (par in names(param)){
      r_loc <- append(r_loc, E(param[[par]]) / post_r_loc[[par]])
      r_scale <- append(r_scale, sd(param[[par]]) / post_r_scale[[par]])
    }
    converged <- all(r_loc > 0.9 && r_loc < 1.1 && r_scale > 0.9 && r_scale < 1.1 )
    return (converged | is.na(converged)) #result NA if param, param_next are same `draws_rvars`
  }else if(type == "rel") {
    return(all(unlist(lapply(target_vars,
           FUN = function(tv) cjs_dist(param[[tv]], param_next[[tv]]) < 0.5 * unlist(dpp_init[tv])))))
  }
}

# Updating parameter value `param` to `param_next` for valid comparison and to control `S`
# Resample `param_next` sample matrix for with PIT weight of F_{prior}F_{post}^{-1}
# S * n_sample posterior S as comparison threshold is possible for the same number of samples
##'
##' @param param numeric vector of prior values i.e. parameter values to be tested
##' @param post_vec posterior samples of length S * M

##' @return resampled posterior with prior information
##' @export
prop_param <-function(param, post_vec){
  S <- niterations(param)
  param_ord <- sort(draws_of(param))
  return (rvar(draws_of(resample_draws(as_draws(rvar(param_ord)),
                                       tabulate(ecdf(param_ord)(post_vec) * S, nbins = S))[[1]])))
}

update_param <-function(param, result, target_vars, cnt = 0, delivDir, thin = thin){
  S <- ndraws(param[[1]])
  M <- ceiling(dim(SBC_fit_to_draws_matrix(result$fits[[1]]))[1] / thin)
  next_param <- param # template
  #next_param_mtr <- matrix(NA, nrow = S, ncol = S) #summarize n_sample to S
  g <- list()
  for (tv in target_vars){
    param_ord <- sort(c(as_draws_df(param)[[tv]]))
    post_mtr <- matrix(NA, nrow = S, ncol = M)
    for (i in 1:S){
      if(!is.null(result$fits[[i]])){ # ADVI can return NULL
        fit_mtr <- subset_draws(SBC_fit_to_draws_matrix(result$fits[[i]]), variable = tv) # change for tv, i order?
        fit_thinned <- posterior::thin_draws(fit_mtr, thin)
        post_mtr[i,] <- c(fit_thinned)
      }
    }
    post_vec <-c(post_mtr)[!is.na(c(post_mtr))]
    next_param_vec <- draws_of(resample_draws(as_draws(rvar(param_ord)), tabulate(ecdf(param_ord)(post_vec) * S, nbins = S))[[1]])
    g[[tv]] <- ppc_hist(param_ord, matrix(next_param_vec, ncol = length(param_ord)))
    next_param[[tv]] <- rvar(next_param_vec)
  }
  ggarrange(g[[target_vars[1]]], g[[target_vars[2]]], nrow = 2)
  ggsave(file =  file.path(delivDir, paste0(paste0(paste0(paste(target_vars[1], target_vars[2]), cnt, "_"), "pp.png", sep = ""))),  bg = "white")
  return (next_param)
}
