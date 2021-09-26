##' Auto calibrate the initial prior samples using SBC iteration
##'
##' @param generator function that generates datasets given each value in `param`
##' @param clamp_dist scalar clamp_disteter values e.g. prior_width
##' @param target_params a named list defining parameter to mixture index mapping
##' @param bandwidth the smoothing bandwidth parameter to pass to stats.density
##' @param mixture_means_init the initial mixture mean rvar
##' @param mixture_bw_init the initial mixture bandwidth rvar
##' @param thin Integer defining thinning parameter
##' @param max_selfcalib_iters the maximum number of iterations to run calibration. if not given will run indefinitely
##' @param nsims_fn function with input: (mixture_means_rvar, mixture_bw_rvar), output: int
##'                 int is future number of parallel datasets to generate given true and its fitted hyperparameter (mixture_means)
##'
##' @param fixed_generator_args *named list* containing additional arguments to generator
##' @export
self_calib <- function(generator, backend, target_params, mixture_means_init, mixture_bw_init, nsims_fn,
                       bandwidth, thin, max_selfcalib_iters, fixed_generator_args){

  if(!posterior::is_rvar(mixture_means_init)){
    mixture_means_init_rvar <- rvar(mixture_means_init)
  }

  if(!posterior::is_rvar(mixture_bw_init)){
    mixture_bw_init_rvar <- rvar(mixture_bw_init)
  }

  if(missing(nsims_fn)){
    nsims_fn <- function(mixture_means_rvar, mixture_means_hat_rvar){10}  # Set default # of SBC iterations to 10
    message(paste("number of simulations has been unspecified, using default value of", nsims))
  }

  if(missing(bandwidth)){
    bandwidth <- "nrd0"
  }

  if(missing(max_selfcalib_iters)){
    max_selfcalib_iters <- Inf
  }

  target_names <- names(target_params) # names of named list, only calibrate these parameters
  ntarget_params <- length(target_names)
  selfcalib_itercount <- 1
  sbc_result <- NULL
  cjs_record <- list()
  for(tp in target_params){
    cjs_record[[tp]] <- c()
  }

  while(selfcalib_itercount < max_selfcalib_iters){
    if(selfcalib_itercount == 1){
      nsims <- 10
      mixture_means_rvar <- mixture_means_init_rvar
      mixture_bw_rvar <- mixture_bw_init_rvar
    }else if(stop){
      message(paste("self_calib terminated on iteration", selfcalib_itercount))
      break
    }else{
      nsims <- nsims_fn(mixture_means_rvar, mixture_means_next_rvar)
      mixture_means_rvar <- mixture_means_next_rvar
      mixture_bw_rvar <- mixture_bw_next_rvar
    }
    dataset <- do.call(generator, c(list(nsims, mixture_means_rvar, mixture_bw_rvar), fixed_generator_args))
    message(paste("Running self-calib iteration", selfcalib_itercount, "with nsims =", nsims))
    sbc_result <- SBC::compute_results(dataset, backend, thin_ranks = thin)
    ndraws <- posterior::ndraws(sbc_result$fits[[1]]$draws())

    mixture_means_hat <- array(rep(NA, ntarget_params * nsims), dim = c(ntarget_params, nsims))
    mixture_bw_hat <- rep(NA, ntarget_params)

    for(tp in target_names){
      pooled_draws <- c()
      for(s in 1:nsims){
        pooled_draws <- c(pooled_draws, posterior::extract_variable(sbc_result$fits[[s]]$draws(), tp))
      }
      gmm_fit <- Mclust(pooled_draws, G = ndataset)
      mixture_means_hat[target_params[[tp]], ] <- gmm_fit$parameters$mean
      mixture_bw_hat[target_params[[tp]], ] <- sqrt(gmm_fit$parameters$variance$sigmasq)
      # fitted_kde <- density(pooled_draws, bw = bandwidth)
      # mixture_means_hat[target_params[[tp]], ] <- pooled_draws[sample(1:(nsims * ndraws), size = nsims)]  # SM2S
      # mixture_bw_hat[target_params[[tp]]] <-fitted_kde$bw
      mixture_means_next <- update_means(mixture_means_hat[target_params[[tp]], ], mixture_means_hat[target_params[[tp]], ])
    }
    mixture_means_hat_rvar <- rvar(array(mixture_means_hat, dim = c(ntarget_params, nsims)))
    mixture_bw_hat_rvar <- rvar(array(mixture_bw_hat, dim = c(nsims, 1)))

    mixture_means_next_rvar <- update_means(mixture_means_rvar, mixture_means_hat_rvar)
    mixture_bw_next_rvar <- update_bw(mixture_bw_rvar, mixture_bw_hat_rvar)

    stop <- TRUE
    for(tp in target_names){
      cjs_record[[tp]] <- c(cjs_record[[tp]], cjs_dist(mixture_means_rvar, mixture_means_hat_rvar[target_params[[tp]], ]))
      if(cjs_record[[tp]][selfcalib_itercount] >= 0.5 * cjs_record[[tp]][1]){
        message(paste("cjs_dist for parameter", tp, ":", cjs_record[[tp]][[selfcalib_itercount]]))
        stop <- FALSE
        break
      }
    }

    selfcalib_itercount <- selfcalib_itercount + 1
  }
  return(sbc_result)
}
# if bw required, fucntion(mixture_means_rvar, mixture_bw_rvar, mixture_mean_hat_rvar, mixture_bw_hat_rvar)
# receive vector for one parameter
update_means <- function(mixture_means_rvar, mixture_means_hat_rvar){
  return(mixture_means_rvar * mean(mixture_means_rvar/ mixture_means_hat_rvar))
}
update_bw <- function(mixture_bw_rvar, mixture_bw_hat_rvar){
  return(mixture_bw_rvar * mean(mixture_bw_rvar/ mixture_bw_hat_rvar))
}

