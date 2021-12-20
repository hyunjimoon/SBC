

##' Auto calibrate the initial prior samples using SBC iterations, with an adaptive update strategy
##'
##' @param generator function that generates datasets given each value in `param`
##' @param backend backend object to use for running SBC
##' @param updator hyperparameter update type
##' @param target_param list of strings indicating target parameter names
##' @param init_mu initial lambda_mu value to use
##' @param init_sigma initial lambda_sigma value to use
##' @param nsims number of datasets i.e. prior draws
##' @param gamma convergence speed e.g. step size
##' @param tol tolerence for determining termination
##' @param fixed_args *named list* containing additional arguments to pass to generator, *after mu and sigma*
##' @export
self_calib_adaptive <- function(generator, backend, updator, target_params, init_lambdas, nsims, gamma, tol, fixed_args){
  dist_types <- fixed_args$dist_types

  ############3
  gamma_estimator <- function(x){
    # mle estimate of shape, scale parameter for gamma distribution
    N <- length(x)
    sum_x <- sum(x)
    sum_log_x <- sum(log(x))
    sum_x_mul_log_x <- sum(x * log(x))

    k_hat <- (N * sum_x) / (N * sum_x_mul_log_x - sum_log_x * sum_x)
    theta_hat <- 1 / N ^ 2 * (N * sum_x_mul_log_x - sum_log_x * sum_x)

    # bias correction
    theta_hat <- N  / (N - 1) * theta_hat
    k_hat <- k_hat - 1 / N * (3 * k_hat - 2/3 * (k_hat / (1 + k_hat)) - 4/5 * (k_hat / (1 + k_hat) ^ 2))

    return(list(alpha=k_hat, beta=1 / theta_hat))
  }

  lognormal_estimator <- function(x){
    n <- length(x)
    mu_hat <- sum(log(x)) / n
    sigma_hat <- sum((log(x) - mu_hat) ^ 2) / (n - 1)
    return(list(mu=mu_hat, sigma=sigma_hat))
  }

  calculate_dap <- function(current_lambdas){
    nsims <- fixed_args$nsims
    datasets <- do.call(generator, list(current_lambdas, fixed_args = fixed_args))
    sbc_result <- SBC::compute_results(datasets, backend, thin_ranks = 1)
    draws_etas <- list()
    return_lambdas <- list()
    for(fit in sbc_result$fits){
      samples <- SBC_fit_to_draws_matrix(fit)
      for(target_param in target_params){

        draws_etas[[target_param]] <- c(draws_etas[[target_param]], posterior::extract_variable(samples, target_param))
      }
    }
    for(target_param in target_params){
      if(fixed_args$dist_type[[target_param]] == "normal"){
        mu <- mean(draws_etas[[target_param]])
        sigma <- sd(draws_etas[[target_param]])
        return_lambdas[[target_param]] <- list(mu=mu, sigma=sigma)
      }
      else if(fixed_args$dist_type[[target_param]] == "gamma"){
        gamma_params <- tryCatch(
          {
            gamma_est = MASS::fitdistr(draws_etas[[target_param]], "gamma")$estimate#, start=list(shape=current_lambdas[[target_param]]$alpha, rate=current_lambdas[[target_param]]$beta))$estimate
            alpha = as.numeric(gamma_est["shape"])
            beta = as.numeric(gamma_est["rate"])
            return(list(alpha=alpha, beta=beta))
          },
          error=function(err){
            message(err)
            message(sprintf("\ngamma mle estimation for parameter %s failed. Falling back to closed form approximation", target_param))
            return(gamma_estimator(draws_etas[[target_param]]))
          }
        )

        return_lambdas[[target_param]] <- gamma_params
      }
      else if(fixed_args$dist_type[[target_param]] == "lognormal"){
        lognormal_params <- tryCatch(
          {
            gamma_est = MASS::fitdistr(draws_etas[[target_param]], "lognormal")$estimate
            mean = as.numeric(gamma_est["meanlog"])
            sd = as.numeric(gamma_est["sdlog"])
            return(list(mu=mean, sigma=sd))
          },
          error=function(err){
            message(err)
            message(sprintf("\nlognormal mle estimation for parameter %s failed. Falling back to closed form approximation", target_param))
            return(lognormal_estimator(draws_etas[[target_param]]))
          }
        )

        return_lambdas[[target_param]] <- lognormal_params
      }
    }
    return(list(return_lambdas = return_lambdas, draws_etas = draws_etas))
  }

  ###############
  # define update strategies

  normal_str_update <- function(draws_dap_lambdas, lambda, gamma){
    normal_str <- function(Tx, x, gamma){
      #b_t <- (Tx + x) / 2
      #Tx + (1/(b_t - Tx) ^2) * (x - Tx)^3
      (Tx + x)/2
    }

    logalpha_new <- log(normal_str(exp(dap$logalpha), exp(lambda$logalpha), gamma))
    logbeta_new <- log(normal_str(exp(dap$logbeta), exp(lambda$logbeta), gamma))

    list(logalpha = logalpha_new, logbeta = logbeta_new)
  }

  mc_update <- function(draws_dap_lambdas, lambdas){
    draws_dap_lambdas
  }

  ###############
  lambda_loss <- function(dap_lambdas, new_lambdas) {
    #return((dap$mu - lambda$mu)^2 + ((dap$sigma)^2 - exp(lambda$logsigma)^2)^2)
    sum((unlist(dap_lambdas) - unlist(new_lambdas))^2)
  }

  eta_loss <- function(dap_eta, new_lambdas) {
    if("mu" %in% names(new_lambdas)){  # normal
      eta <- rnorm(length(dap_eta), mean=new_lambdas$mu, sd = new_lambdas$sigma)
    }
    else if("alpha" %in% names(new_lambdas)){  # gamma
      eta <- rgamma(length(dap_eta), shape=new_lambdas$alpha, rate = new_lambdas$beta)
    }
    return(cjs_dist(eta, dap_eta))
  }

  normal_kl_divergence <- function(dap, lambda){
    v_1 <- dap$sigma^2
    v_2 <- lambda$sigma^2
    (dap$mu - lambda$mu)^ 2/(2 * v_2) + 0.5 * ((v_1 / v_2) - log(v_1 / v_2) - 1)
  }

  # end function declarations
  lambda_current <- init_lambdas
  t_df <- list()
  for (iter_num in 1:niter) {
    stop <- TRUE
    t_df$iter <- c(t_df$iter, iter_num)

    dap_result <- calculate_dap(lambda_current)

    dap_tx_plot_list <- list()
    dap_tx_plot_index = 1
    for(target_param in target_params){
      param_lambdas <- lambda_current[[target_param]]
      lambda_count <- length(param_lambdas)
      for(i in 1:lambda_count){
        lambda_colname <- paste(target_param, names(param_lambdas)[i], sep = "_")
        t_df[[lambda_colname]] <- c(t_df[[lambda_colname]], param_lambdas[[names(param_lambdas)[i]]])
      }
      if (dist_types[[target_param]] == "normal"){
        prior_dist_samples <- rnorm(length(dap_result$draws_etas[[target_param]]), mean=param_lambdas$mu, sd=param_lambdas$sigma)
      }
      else if (dist_types[[target_param]] == "gamma"){
        prior_dist_samples <- rgamma(length(dap_result$draws_etas[[target_param]]), shape=param_lambdas$alpha, r=param_lambdas$beta)
      }
      else if (dist_types[[target_param]] == "lognormal"){
        prior_dist_samples <- rlnorm(length(dap_result$draws_etas[[target_param]]), meanlog = param_lambdas$mu, sdlog = param_lambdas$sigma)
      }

      plot_df <- data.frame(dap=dap_result$draws_etas[[target_param]], prior=prior_dist_samples)
      plot <- ggplot2::ggplot(plot_df) + ggplot2::geom_density(aes(x=dap), color="red") + ggplot2::geom_density(aes(x=prior)) + ggplot2::ggtitle(sprintf("%s (red=dap)", target_param))
      if(iter_num == 1 || iter_num %% 10 == 0){
        print(plot)
      }
      dap_tx_plot_list[[dap_tx_plot_index]] <- plot
      dap_tx_plot_index <- dap_tx_plot_index + 1
    }

    if(updator == "normal_str_update"){
      stop("Unfinished implementation")
      lambda_new <- normal_str_update(dap_result$return_lambdas, lambda_current, gamma)
    }
    else if(updator == "mc_update"){
      lambda_new <- mc_update(dap_result$return_lambdas, lambda_current)
    }

    message(sprintf("Iteration %d:", iter_num))
    for(target_param in target_params){
      param_lambdas <- lambda_current[[target_param]]

      param_lambda_loss <- lambda_loss(dap_result$return_lambdas[[target_param]], param_lambdas)
      param_eta_loss <-  eta_loss(dap_result$draws_etas[[target_param]], param_lambdas)

      t_df[[paste(target_param, "lambda_loss", sep="_")]] <- c(t_df[[paste(target_param, "lambda_loss", sep="_")]], param_lambda_loss)
      t_df[[paste(target_param, "eta_loss", sep="_")]] <- c(t_df[[paste(target_param, "eta_loss", sep="_")]], param_eta_loss)
      message(sprintf("parameter %s - lambda loss: %f eta_loss: %f", target_param, param_lambda_loss, param_eta_loss))
      if(all(abs(unlist(param_lambdas) - unlist(dap_result$return_lambdas[[target_param]])) < tol) && iter_num > 1){
        stop <- TRUE && stop
      }
      else{
        stop <- FALSE
      }
    }

    if(stop){
      message(sprintf("Terminating self_calib on iteration %d", iter_num))
      break
    }
    lambda_current <- lambda_new
  }
  t_df <- as.data.frame(t_df)
  return(list(lambda=lambda_current, t_df=t_df))
}



##' Auto calibrate the initial prior samples using SBC iteration for gaussian approximation
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
self_calib_gaussian <- function(generator, backend, mixture_means_init_draws_rvars, mixture_bw_init_draws_rvars, nsims_fn,
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
    datasets <- do.call(generator, c(list(mixture_means_draws_rvars, mixture_bw_draws_rvars), fixed_generator_args))
    message("generator returned value")
    sbc_result <- SBC::compute_results(datasets, backend, thin_ranks = thin)
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
      mixture_means_next_draws_rvars[[target_param_name]] <- update_quantile_approximation(mixture_means_draws_rvars[[target_param_name]], posterior::resample_draws(posterior::rvar(array(rep(as.vector(gmm_fit$parameters$mean), each = nsims), dim = c(nsims, nsims))), weights = prop_est))
      mixture_bw_next_draws_rvars[[target_param_name]] <- update_quantile_approximation(mixture_means_next_draws_rvars[[target_param_name]])
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

##' Auto calibrate the initial prior samples using SBC iteration and gmm approximation
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
  self_calib_gmm <- function(generator, backend, mixture_means_init_draws_rvars, mixture_bw_init_draws_rvars, nsims_fn,
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
    datasets <- do.call(generator, c(list(mixture_means_draws_rvars, mixture_bw_draws_rvars), fixed_generator_args))
    message("generator returned value")
    sbc_result <- SBC::compute_results(datasets, backend, thin_ranks = thin)
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
      mixture_means_next_draws_rvars[[target_param_name]] <- update_quantile_approximation(mixture_means_draws_rvars[[target_param_name]], posterior::resample_draws(posterior::rvar(array(rep(as.vector(gmm_fit$parameters$mean), each = nsims), dim = c(nsims, nsims))), weights = prop_est))
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
#' @return vector of phi which are quantile function values
approx_quantile_phi <- function(draws, S) {
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

#' Update mixture mean through quantile approxiation based on analytic gradient of quantile loss
#'
#' @param hyperparam_rvar a posterior::rvar object of prior(pre-transformation) mixture mean values
#' @param hyperparam_hat_rvar a posterior::rvar object of posterior(post-transformation) mixture mean values
#' @param S Number of approximation points for the target quantile function.
#' @param n_post_samples the number of samples to draw from posterior, to approximate the expected quantile loss
#' @param epsilon gradient update coefficient
#' @return a posterior::rvar object with the same dimension as the input rvars.
#' @export
update_quantile_approximation <- function(hyperparam_rvar, hyperparam_hat_rvar, S,  n_post_samples,  epsilon) {
  phi <- approx_quantile_phi(hyperparam_rvar, S = S)
  phi_post <- approx_quantile_phi(hyperparam_hat_rvar, S = S)
  updated_phi <- phi
  wass_last <- wasserstein(updated_phi, phi_post)
  iters <- 0
  while(iters <= 5 || abs(wasserstein(updated_phi, phi_post) - wass_last) > wass_last * 0.001){
    wass_last <- wasserstein(updated_phi, phi_post)
    #plot(phi, unlist(lapply(c(1:S), function(x) {(2 * x - 1) / (2 * S)})), type = "l")
    #lines(updated_phi, unlist(lapply(c(1:S), function(x) {(2 * x - 1) / (2 * S)})), col="red")
    zprime <- sample_quantile_phi(n_post_samples, updated_phi)
    for(s in 1:S) {
      # delta_sum <- 0
      # for(m in 1:n_post_samples){
      #   zprime_delta <- sum(zprime[m] < updated_phi[s])
      #   delta_sum <- delta_sum + ((2 * s - 1) / (2 * S) - zprime_delta / n_post_samples + if (zprime[m] - updated_phi[s] == 0) 1 else 0)
      # }
      delta_sum <-(2 * s - 1) / (2 * S) - (sum(zprime < updated_phi[s]) + sum(zprime == updated_phi[s])) / n_post_samples
      updated_phi[s] <- updated_phi[s] + epsilon * (delta_sum / n_post_samples)
      #zprime <- sample_quantile_phi(n_post_samples, updated_phi)
      #zprime_delta <- sum(zprime < updated_phi[s])
      #print(paste(zprime_delta / n_post_samples, "/", (2 * s - 1) / (2 * S)))
      #updated_phi[s] <- updated_phi[s] + epsilon * ((2 * s - 1) / (2 * S) - zprime_delta / n_post_samples)  # (tau_{i - 1} + tau_i) / S = (s / S + (s - 1) / S) / 2

      #phi_delta <- c(phi_delta, ((2 * s - 1) / (2 * S) - zprime_delta / n_post_samples))
    }
    iters <- iters + 1
    #message(updated_phi)
    #message(paste(wasserstein(updated_phi, phi_post), wass_last))
  }
  print(paste("optimization iters:", iters))
  return(list(updated_phi=posterior::rvar(array(rep(updated_phi, each = nsims), dim = c(nsims, nsims))), phi=phi)) # currently all nsims receive same updated mus
}
