##' Auto calibrate the initial prior samples using SBC iteration
##'
##' @param generator function that generates datasets given each value in `param`
##' @param clamp_dist scalar clamp_disteter values e.g. prior_width
##' @param target_variables a named list defining parameter to mixture index mapping
##' @param bandwidth the smoothing bandwidth parameter to pass to stats.density
##' @param mixture_mean_init the initial mixture mean rvar
##' @param mixture_bw_init the initial mixture bandwidth rvar
##' @param thin Integer defining thinning parameter
##' @param max_calib_iterations the maximum number of iterations to run calibration. if not given will run indefinitely
##' @param sbc_iterations An integer returning the number of SBC iterations or a function with the signature (mixture_mean_rvar, mixture_bw_rvar) -> int
##'                       which returns the number of SBC iterations to run given previous calibration run results.
##' @param fixed_generator_args *named list* containing additional arguments to generator
##' @export
self_calib <- function(generator, backend, target_variables, bandwidth, mixture_mean_init, mixture_bw_init, thin, max_calib_iterations, sbc_iterations, fixed_generator_args){
  if(missing(sbc_iterations)){
    sbc_iterations <- 10  # Set default # of SBC iterations to 10
    message(paste("sbc_iterations has been unspecified, using default value of", sbc_iterations))
  }
  if(missing(max_calib_iterations)){
    max_calib_iterations <- Inf
  }
  if(missing(bandwidth)){
    bandwidth <- "nrd0"
  }
  target_var_names <- names(target_variables) # names of named list, only run calibration for the following parameters
  target_var_count <- length(target_var_names)

  if(!posterior::is_rvar(mixture_mean_init)){
    mixture_mean_init <- rvar(mixture_mean_init)
  }

  initial_sbc_iterations <- sbc_iterations(mixture_mean_init, mixture_bw_init)
  initial_mixture_mean_rvar <- mixture_mean_init

  draw_from_mixture <- function(mm_rvar, n_samples, n_sbc_iters){
    # draw from a rvar of dim (M, N), draw = 1 along row(M), each n_samples
    mm_mat <- draws_of(mm_rvar)[1, , ]
    ret_vec <- c()
    S <- n_sbc_iters
    N <- dim(mm_mat)[1]
    for(m in 1:n_samples){  # I truly, truly hate R array indexing
      for(n in 1:N){
        for(S in 1:S){

          ret_vec <- c(ret_vec, sample(mm_mat[n, ], 1, replace = TRUE))
        }
      }
    }
    rvar(array(ret_vec, dim = c(n_sbc_iters, N, n_samples)))
  }
  #mixture_mean_rvar <- draw_from_mixture(mixture_mean_init, dim(mixture_mean_init)[2], initial_sbc_iterations)
  dataset <- do.call(generator, c(list(initial_sbc_iterations, initial_mixture_mean_rvar, mixture_bw_init), fixed_generator_args))

  iteration_counter <- 1
  previous_draws <- NULL
  sbc_result <- NULL

  cjs_record <- list()
  for(tv in target_variables){
    cjs_record[[tv]] <- c()
  }
  while(iteration_counter < max_calib_iterations){

    S <- dim(dataset$parameters)[[1]]  # number of SBC iterations
    message(paste("Running self_calib iteration", iteration_counter, "with S =", S))
    sbc_result <- SBC::compute_results(dataset, backend, thin_ranks = thin)
    n_posterior_draws <- posterior::ndraws(sbc_result$fits[[1]]$draws())


    extracted_mixture_means <- array(rep(NA, target_var_count * S * n_posterior_draws), dim = c(target_var_count, S * n_posterior_draws))
    extracted_mixture_bw <- rep(NA, target_var_count)
    for(tv in target_var_names){
      pooled_draws <- c()
      for(s in 1:S){
        pooled_draws <- c(pooled_draws, posterior::extract_variable(sbc_result$fits[[s]]$draws(), tv))
      }
      fitted_kde <- density(pooled_draws, bw = bandwidth)
      #theta_tf_next[[tv]] <- rvar(array(sample(kde_fit_points, S, replace = TRUE), dim = c(S, 1))) # draw from the posterior distribution
      extracted_mixture_bw[target_variables[[tv]]] <- fitted_kde$bw
      extracted_mixture_means[target_variables[[tv]], ] <- pooled_draws
      cjs_record[[tv]] <- c(cjs_record[[tv]], cjs_dist(as.vector(posterior::draws_of(initial_mixture_mean_rvar[target_variables[[tv]], ])), pooled_draws))
    }

    extracted_mixture_mean_rvar <- rvar(array(extracted_mixture_means, dim = c(S, target_var_count, S * n_posterior_draws)))
    extracted_mixture_bw_rvar <- rvar(array(extracted_mixture_bw, dim = c(S, 1)))

    stop <- TRUE
    for(tv in target_var_names){
      if(cjs_record[[tv]][iteration_counter] >= 0.5 * cjs_record[[tv]][1]){
        message(paste("cjs_dist for parameter", tv, ":", cjs_record[[tv]][[iteration_counter]]))
        stop <- FALSE
        break
      }
    }
    if(stop && iteration_counter != 1){
      message(paste("self_calib terminated on iteration", iteration_counter))
      break
    }
    previous_draws <- as.vector(extracted_mixture_means)


    dataset <- do.call(generator, c(list(sbc_iterations(mixture_mean_init, mixture_bw_init), extracted_mixture_mean_rvar, extracted_mixture_bw_rvar),
                                    fixed_generator_args))
    iteration_counter <- iteration_counter + 1
  }
  return(sbc_result)
}


