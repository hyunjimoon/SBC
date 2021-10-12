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
      prop_est <- gmm_fit$parameters$pro
      # vectorize gmmmeans -> sample ->
      mixture_means_next_draws_rvars[[target_param_name]] <- update_means(mixture_means_draws_rvars[[target_param_name]], posterior::resample_draws(posterior::rvar(array(rep(as.vector(gmm_fit$parameters$mean), each = nsims), dim = c(nsims, nsims))), weights = prop_est))
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
################
# quantile approximation

#' Given a vector of draws, return a vector of length S of phi which best approximates the CDF.
#' @param draws a vector of sample draws from a distribution
#' @param S number of quantile points
#' @return vector of phi where phi in \[0, 1\]
approx_quantile_phi <- function(draws, S) {
  probs <- c(1:S)
  probs <- unlist(lapply(c(1:S), function(x) {(2 * x - 1) / (2 * S)}))  # generate (tau_i + tau_{i+1})/2
  return(quantile(draws, probs, names = FALSE))
}

#' Given a vector of phis, which represent the quantiles of equally spaced probabilities on [0, 1] defined as i/S, return a function that returns N random samples from the quantile function
#' @param N number of samples to draw, returned as a vector of length N
#' @param phis vector of phis to sample from
#' @return a vector of samples drawn by inverse transform sampling
sample_quantile_phi <- function(N, phis) {
  return(sample(phis, N, replace = TRUE))
}


#' Calculate qualtile huber loss for a given tau tau_{s_index}
#' @param phi_prior vector of phis for the prior(pre-transformation) distribution
#' @param phi_post vector of phis for the posterior(post-transformation) distribution
#' @param s_index The tau index to calculate loss
#' @param k interval to calculate loss. resulting loss will be in range \[-k, k\]
#' @param S the number of phis. equal to length(phi_prior) = length(phi_post)
#' @param n_post_samples the number of samples to draw from posterior, to approximate the expected huber loss
#' @return a vector of length 2, where the first value is the expected huber loss and the second value the number of posterior samples less than the quantile value at phi\[s_index\]
quantile_huber_loss <- function(phi_prior, phi_post, s_index, k, S, n_post_samples) {
  summed_rho_mean <- 0
  zprime_delta <- 0
  for(n in 1:n_post_samples){
    zprime <- sample(phi_post, 1)
    u <- zprime - phi_prior[s_index]
    delta <- sum(u < 0)
    tau_hat <- s_index / S
    huber_loss <- if (abs(u) <= k ) 1/2 * u ** 2  else k * (abs(u) - 1/2 * k)
    rho <- abs(tau_hat - u) * huber_loss / k

    summed_rho_mean <- summed_rho_mean + rho
    zprime_delta <- zprime_delta + if(zprime < phi_prior[s_index]) 1 else 0
  }
  return(c(summed_rho_mean / n_samples, zprime_delta))
}


#' Update mixture mean through quantile approxiation based on huber loss
#'
#' @param mixture_means_rvar a posterior::rvar object of prior(pre-transformation) mixture mean values
#' @param mixture_means_rvar a posterior::rvar object of posterior(post-transformation) mixture mean values
#' @param S Number of approximation points for the target quantile function.
#' @param k the huber loss bound. huber loss returns in the boundary of \[-k, k\]
#' @param n_post_samples the number of samples to draw from posterior, to approximate the expected huber loss
#' @param epsilon gradient update coefficient
#' @return a posterior::rvar object with the same dimension as the input rvars.
#' @export
update_quantile_approximation <- function(mixture_means_rvar, mixture_means_hat_rvar, S, k, n_post_samples,  epsilon) {
  phi <- approx_quantile_phi(posterior::draws_of(mixture_means_rvar), S = S)
  phi_post <- approx_quantile_phi(posterior::draws_of(mixture_means_hat_rvar), S = S)
  updated_phi <- phi
  summed_loss <- 0
  for(s in 1:S) {
    qhl = quantile_huber_loss(phi, phi_post, s, k, S, n_post_samples)
    huber_loss <- qhl[1]
    zprime_delta <- qhl[2]
    summed_loss <- summed_loss + huber_loss
    updated_phi[s] <- updated_phi[s] + epsilon * (s / S - zprime_delta)
  }
  posterior::rvar(array(rep(sample_quantile_phi(nsims, updated_phi), each = nsims), dim = c(nsims, nsims)))  # currently all nsims receive same updated mus
}

