##' Auto calibrate the initial prior samples using SBC iteration
##'
##' @param generator function that generates datasets given each value in `param`
##' @param backend A backend object to use for running SBC
##' @param mixture_means_init_draws_rvars the initial mixture mean draws_rvars
##' @param mixture_bw_init_draws_rvars the initial mixture bandwidth draws_rvars
##' @param nsims_fn function with input: (mixture_means_rvar, mixture_bw_rvar), output: int
##'                 int is future number of parallel datasets to generate given true and its fitted hyperparameter (mixture_means)
##' @param thin Integer defining thinning parameter
##' @param max_selfcalib_iters the maximum number of iterations to run calibration. if not given will run indefinitely
##' @param save_all_results Boolean if TRUE returns a list of all SBC results, FALSE returns just the result of the last iteration.
##' @param transform_types Transformtype for mixture fitting
##' @param fixed_generator_args *named list* containing additional arguments to pass to generator, *after mixture_means_draws_rvars and mixture_bw_draws_rvars*
##' @export
self_calib <- function(generator, backend, mixture_means_init_draws_rvars, mixture_bw_init_draws_rvars, nsims_fn,
                        thin, max_selfcalib_iters, save_all_results, transform_types, fixed_generator_args){
  if(missing(nsims_fn)){
    nsims_fn <- function(...){300}  # Set default # of SBC iterations to 30
    message(paste("number of simulations has been unspecified, using default value of", nsims))
  }
  target_params <- posterior::variables(mixture_means_init_draws_rvars) # names of named list, only run calibration for the following parameters
  ntarget_params <- length(target_params)

  if(missing(max_selfcalib_iters)){
    max_selfcalib_iters <- Inf
  }

  if(missing(save_all_results)){
    save_all_results <- FALSE
  }

  selfcalib_itercount <- 1
  sbc_result <- NULL
  cjs_record <- list()
  for(tp in target_params){
    cjs_record[[tp]] <- c()
  }

  if(save_all_results){
    sbc_result_env <- new.env()
  }

  while(selfcalib_itercount < max_selfcalib_iters){
    if(selfcalib_itercount == 1){
      nsims <- nsims_fn(1)
      mixture_means_draws_rvars <- mixture_means_init_draws_rvars
      mixture_bw_draws_rvars <- mixture_bw_init_draws_rvars
    }else{
      nsims <- nsims_fn(mixture_means_draws_rvars, mixture_means_next_draws_rvars)
      mixture_means_draws_rvars <- mixture_means_next_draws_rvars
      mixture_bw_draws_rvars <- mixture_bw_next_draws_rvars
    }
    message(paste("Running self-calib iteration", selfcalib_itercount, "with nsims =", nsims))
    message("Calling generator..")
    dataset <- do.call(generator, c(list(mixture_means_draws_rvars, mixture_bw_draws_rvars), fixed_generator_args))
    message("generator returned value")
    sbc_result <- SBC::compute_results(dataset, backend, thin_ranks = thin)
    if(save_all_results){
      sbc_result_env[[paste0("result_", selfcalib_itercount)]] <- sbc_result
    }
    for(s in 1:nsims){
      returned_errors <- sbc_result$errors
      if(!is.null(returned_errors[[s]])){
        message("SBC returned 1 or more errors. Terminating and returning the last unsuccessful SBC result...")
        return(sbc_result)
      }
    }
    ndraws <- posterior::ndraws(sbc_result$fits[[1]]$draws())

    mixture_means_next_draws_rvars <- list()
    mixture_bw_next_draws_rvars <- list()
    for(n_variable in 1:ntarget_params){
      target_param_name <- target_params[[n_variable]]
      pooled_draws <- c()
      for(s in 1:nsims){
        pooled_draws <- c(pooled_draws, posterior::extract_variable(sbc_result$fits[[s]]$draws(), target_param_name))
      }
      #pooled_draws = tf_param_vec(pooled_draws, if (missing(transform_types)) "identity" else transform_types[[target_param_name]])
      pooled_draws = tf_param_vec(pooled_draws, transform_types[[target_param_name]])

      gmm_fit <- mclust::Mclust(pooled_draws, G = nsims, verbose = FALSE)
      # vectorize gmmmeans -> sample ->
      mixture_means_next_draws_rvars[[target_param_name]] <- update_means(mixture_means_draws_rvars[[target_param_name]], posterior::rvar(array(rep(as.vector(gmm_fit$parameters$mean), each = nsims), dim = c(nsims, nsims))))
      mixture_bw_next_draws_rvars[[target_param_name]] <- update_bw(mixture_means_next_draws_rvars[[target_param_name]])
    }
    mixture_means_next_draws_rvars <- do.call(draws_rvars, mixture_means_next_draws_rvars)
    mixture_bw_next_draws_rvars <- do.call(draws_rvars, mixture_bw_next_draws_rvars)

    stop <- TRUE
    for(tp in target_params){
      cjs_record[[tp]] <- c(cjs_record[[tp]], cjs_dist(mixture_means_draws_rvars[[tp]], mixture_means_next_draws_rvars[[tp]]))
      if(cjs_record[[tp]][selfcalib_itercount] >= 0.5 * cjs_record[[tp]][1]){
        message(paste("cjs_dist for parameter", tp, ":", cjs_record[[tp]][[selfcalib_itercount]]))
        stop <- FALSE
      }
    }
    if(stop){
      message(paste("self_calib terminated on iteration", selfcalib_itercount))
      break
    }
    selfcalib_itercount <- selfcalib_itercount + 1
  }
  return(if(save_all_results) sbc_result_env else sbc_result)
}
# fucntion(mixture_means_rvar, mixture_bw_rvar, mixture_mean_hat_rvar, mixture_bw_hat_rvar) possible
update_means <- function(mixture_means_rvar, mixture_means_hat_rvar){
  return(mixture_means_rvar * mean(mixture_means_rvar/ mixture_means_hat_rvar))
}
update_bw <- function(mixture_means_next_rvars){
  return (rvar(rep(bw.nrd0(draws_of(mixture_means_next_rvars)), posterior::niterations(mixture_means_next_rvars))))
}
