

##' Auto calibrate the initial prior samples using SBC iterations, with an adaptive update strategy
##'
##' @param generator function that generates datasets given each value in `param`
##' @param backend backend object to use for running SBC
##' @param updator hyperparameter update type
##' @param target_param string type target parameter name
##' @param init_mu initial lambda_mu value to use
##' @param init_sigma initial lambda_sigma value to use
##' @param nsims number of datasets i.e. prior draws
##' @param gamma convergence speed e.g. step size
##' @param tol tolerence for determining termination
##' @param fixed_args *named list* containing additional arguments to pass to generator, *after mu and sigma*
##' @export
self_calib_adaptive <- function(generator, backend, updator, target_param, init_mu, init_sigma, nsims, gamma, tol, fixed_args){

  calculate_dap <- function(mu, sigma){
    nsims <- fixed_args$nsims
    datasets <- do.call(generator, c(list(mu, sigma), list(fixed_args = fixed_args)))
    sbc_result <- SBC::compute_results(datasets, backend, thin_ranks = 1)
    draws_eta <- c()
    for(fit in sbc_result$fits){
      samples <- SBC_fit_to_draws_matrix(fit)
      draws_eta <- c(draws_eta, posterior::extract_variable(samples, target_param))
    }
    if(is.element("rvar", class(init_mu))){
      gmm_fit <- mclust::Mclust(draws_eta, G = nsims, verbose = FALSE)
      prop_est <- gmm_fit$parameters$pro
      mu <- posterior::rvar(array(rep(sample(as.vector(gmm_fit$parameters$mean), prob = prop_est), each = nsims), dim = c(nsims, nsims)))
      sigma <- rvar(rep(bw.nrd0(draws_of(mu)), posterior::niterations(mu)))
    }else{
      # assume normal for dap
      mu <- mean(draws_eta)
      sigma <- sd(draws_eta)
    }
    return(list(mu=mu, sigma=sigma, draws_eta=draws_eta)) # draws_eta
  }

  ###############
  # define update strategies

  heuristic_update_cubic <- function(dap, lambda, max_diff_mu, max_diff_sigma){
    logsigma_Txgx <- function(Tx, x) {
      (x^2 + x^2) / (Tx + x)
    }
    logsigma_xgTx <- function(Tx, x) {
      (Tx^2 + x^2) / (2*Tx)
    }
    if(dap$logsigma > lambda$logsigma){
      print("T_logsigma > logsigma")
      logsigma_new <- logsigma_Txgx(dap$logsigma, lambda$logsigma)
    }
    else{
      print("T_logsigma <= logsigma")
      logsigma_new <- logsigma_xgTx(dap$logsigma, lambda$logsigma)
    }

    cubic <- function(Tx, x, b_t){
      Tx + 1/(Tx - b_t) ^2 * (x - Tx)^3
      #1/2 * (x - Tx) + Tx
    }

    mu_new <- cubic(dap$mu, lambda$mu, dap$mu + max_diff_mu)
    logsigma_new <- cubic(dap$logsigma, lambda$logsigma, dap$logsigma + max_diff_sigma)
    print(sprintf("T_x:%f x:%f mu_new:%f b_t:%f", dap$mu, lambda$mu, mu_new, dap$mu + max_diff_mu))

    list(mu = mu_new, logsigma = logsigma_new)
  }

  normal_str_update <- function(dap, lambda, gamma){
    normal_str <- function(Tx, x, gamma){
      b_t <- (gamma +1) * Tx + gamma * x
      Tx + (1/(b_t - Tx) ^2) * (x - Tx)^3
    }
    mu_new <- normal_str(dap$mu, lambda$mu, gamma)
    sigmasq_new <- normal_str(exp(dap$logsigma)^2, exp(lambda$logsigma)^2, gamma)
    logsigma_new <- log(sqrt(sigmasq_new))
    print(sprintf("T_x:%f x:%f mu_new:%f", dap$mu, lambda$mu, mu_new))
    print(sprintf("T_x:%f x:%f sigma_new:%f", exp(dap$logsigma)^2, exp(lambda$logsigma)^2, sigmasq_new))

    list(mu = mu_new, logsigma = logsigma_new)
  }

  heuristic_update <- function(dap, lambda){
    logsigma_Txgx <- function(Tx, x) {
      (x^2 + x^2) / (Tx + x)
    }
    logsigma_xgTx <- function(Tx, x) {
      (Tx^2 + x^2) / (2*Tx)
    }
    mu_abs_Txgx <- function(Tx, x) {
      alpha = 1 / sqrt(2) * abs(Tx - x)
      (2 * (abs(x) + abs(alpha)) ^ 2) / (abs(Tx) * abs(x) + 2 * abs(alpha))
      #(x^2 + x^2) / abs(Tx + x)
    }
    mu_abs_xgTx <- function(Tx, x) {
      #(Tx^2 + x^2) / abs(2*Tx)
      alpha = 1 / sqrt(2) * abs(Tx - x)
      (2 * (abs(x) + abs(alpha)) ^ 2) / (abs(Tx) * abs(x) + 2 * abs(alpha))
    }
    if(dap$logsigma > lambda$logsigma){
      print("T_logsigma > logsigma")
      logsigma_new <- logsigma_Txgx(dap$logsigma, lambda$logsigma)
    }
    else{
      print("T_logsigma <= logsigma")
      logsigma_new <- logsigma_xgTx(dap$logsigma, lambda$logsigma)
    }

    if(abs(dap$mu) > abs(lambda$mu)){
      print("T_mu > mu")
      mu_new <- mu_abs_Txgx(dap$mu, lambda$mu) * sign(lambda$mu)
    }
    else{
      print("T_mu <= mu")
      mu_new <- mu_abs_xgTx(dap$mu, lambda$mu) * sign(lambda$mu)
    }

    draws_eta <- dap$draws_eta
    hist(draws_eta, breaks=80, freq = FALSE)
    xval <- seq(min(draws_eta), max(draws_eta), length.out = 100)
    lines(xval, dnorm(xval, dap$mu, dap$sigma))
    lines(xval, dnorm(xval, mu_new, exp(logsigma_new)), col="red")
    print(sprintf("Tx: %f x: %f, new logsigma: %f", dap$logsigma, lambda$logsigma, logsigma_new))
    print(sprintf("Tx: %f x: %f, new mu: %f", dap$mu, lambda$mu, mu_new))
    #print(paste(mu_new, logsigma_new))

    list(mu = mu_new, logsigma = logsigma_new)
  }

  gradient_update <- function(dap, lambda){
    loss <- function(lambda) {
      dap_lambda <- c(dap$mu, dap$sigma)
      message("loss : ", sum((lambda - dap_lambda)^2))
      sum((lambda - dap_lambda)^2)
    }
    grad_loss <- function(lambda) {
      # rough finite difference gradients such that steps are bigger
      # than the expected error caused by the simulations from prior and posterior
      numDeriv::grad(loss, lambda, method.args=list(eps=0.3, d = 0.3))
    }
    # gradient descent update
    # gamma <- 0.5
    lambda <- lambda - gamma * grad_loss(lambda)
    return(list(mu=lambda$mu, logsigma = lambda$logsigma))
  }
  ########
  # mixture update strategy

  quantile_update <- function(dap, lambda){
    phi <- approx_quantile_phi(posterior::draws_of(lambda$mu), S = fixed_args$nsims)
    updated_phi <- phi
    phi_post <- approx_quantile_phi(posterior::draws_of(dap$mu), S = fixed_args$nsims)
    wass_last <- wasserstein(updated_phi, phi_post)
    iters <- 0
    while(iters <= 5 || abs(wasserstein(updated_phi, phi_post) - wass_last) > wass_last * 0.001){
      wass_last <- wasserstein(updated_phi, phi_post)
      #plot(phi, unlist(lapply(c(1:S), function(x) {(2 * x - 1) / (2 * S)})), type = "l")
      #lines(updated_phi, unlist(lapply(c(1:S), function(x) {(2 * x - 1) / (2 * S)})), col="red")
      zprime <- sample_quantile_phi(fixed_args$ndraws, updated_phi)
      S <- fixed_args$nsims
      for(s in 1:S) {
        delta_sum <-(2 * s - 1) / (2 * S) - (sum(zprime < updated_phi[s]) + sum(zprime == updated_phi[s])) / fixed_args$ndraws
        updated_phi[s] <- updated_phi[s] + 0.1 * (delta_sum / fixed_args$ndraws)
      }
      iters <- iters + 1
      #message(updated_phi)
      message(paste(wasserstein(updated_phi, phi_post), wass_last))
    }
    #message(paste("optimization iters:", iters))
    return_mu <- posterior::rvar(array(rep(updated_phi, each = nsims), dim = c(nsims, nsims)))

    return_logsigma <- log(rvar(rep(bw.nrd0(draws_of(return_mu)), posterior::niterations(return_mu))))
    return(list(mu=return_mu, logsigma = return_logsigma)) # currently all nsims receive same updated mus
  }

  ###############
  lambda_loss <- function(dap, lambda) {
    return((dap$mu - lambda$mu)^2 + ((dap$sigma)^2 - exp(lambda$logsigma)^2)^2)
  }

  eta_loss <- function(dap_eta, lambda) {
    eta <- rnorm(length(dap_eta), lambda$mu, exp(lambda$logsigma))
    return(cjs_dist(eta, dap_eta))
  }

  normal_kl_divergence <- function(dap, lambda){
    v_1 <- sqrt(dap$sigma)
    v_2 <- sqrt(lambda$sigma)
    (dap$mu - lambda$mu)^ 2/(2 * v_2) + 0.5 * ((v_1 / v_2) - log(v_1 / v_2) - 1)
  }

  # end function declarations

  mu_current <- init_mu
  sigma_current <- init_sigma
  cjs_prev <- Inf
  t_df <- list(iter=c(), T_logsigma=c(), logsigma=c(), new_logsigma=c(), T_mu=c(), mu=c(), new_mu=c(), lambda_loss=c(), eta_loss=c())
  heuristic_max_diff_mu <- -Inf
  heuristic_max_diff_logsigma <- -Inf
  for (iter_num in 1:niter) {
    stop <- FALSE
    t_df$iter <- c(t_df$iter, iter_num)

    dap_result <- calculate_dap(mu_current, sigma_current)
    dap_result$logsigma = log(dap_result$sigma)

    t_df$mu <- c(t_df$mu, mu_current)
    t_df$logsigma <- c(t_df$logsigma, log(sigma_current))
    t_df$T_mu <- c(t_df$T_mu, dap_result$mu)
    t_df$T_logsigma <- c(t_df$T_logsigma, dap_result$logsigma)

    if(is.element("rvar", class(init_mu))){
      lambda_current <- list(mu=mu_current, logsigma=log(sigma_current))
      mixture_means_next_draws_rvars <- quantile_update(dap_result, lambda_current)
      mu_new <- mixture_means_next_draws_rvars$mu
      logsigma_new <- mixture_means_next_draws_rvars$logsigma
      cjs_record <- cjs_dist(mu_current, mu_new)
      if(iter_num ==1){
        cjs_prev <- cjs_record
      }
      message(sprintf("Iteration %d - cjs_dist for parameter %f", iter_num,cjs_record))
      if(cjs_record >= 0.5 * cjs_prev){
        stop <- TRUE
      }
    }else{
      lambda_current <- list(mu=mu_current, logsigma=log(sigma_current), sigma=sigma_current)

      plot_df <- data.frame(dap=dap_result$draws_eta, prior=rnorm(length(dap_result$draws_eta), mu_current, sigma_current))
      mx <- ggplot(plot_df)
      mx <- mx  + geom_density(aes(x=dap), color="red") + geom_density(aes(x=prior)) + ggtitle("red = dap")
      print(mx)
      ggsave(sprintf("iter_%d.png", iter_num))

      if(updator == "gradient"){
        lambda_new <- gradient_update(dap_result, lambda_current)
      }else if(updator == "heuristic"){
        lambda_new <- heuristic_update(dap_result, lambda_current)
      }
      else if(updator == "heuristic_cubic"){
        heuristic_max_diff_mu <- max(heuristic_max_diff_mu, abs(dap_result$mu - lambda_current$mu) * 1.1)
        heuristic_max_diff_logsigma <- max(heuristic_max_diff_logsigma, abs(dap_result$logsigma - lambda_current$logsigma) * 1.1)
        lambda_new <- heuristic_update_cubic(dap_result, lambda_current, heuristic_max_diff_mu, heuristic_max_diff_logsigma)
      }else if(updator == "normal_str_update"){
        lambda_new <- normal_str_update(dap_result, lambda_current, gamma)
      }

      t_df$new_mu <- c(t_df$new_mu, lambda_new$mu)
      t_df$new_logsigma <- c(t_df$new_logsigma, lambda_new$logsigma)

      mu_new <- lambda_new$mu
      sigma_new <- exp(lambda_new$logsigma)
      if (abs(mu_current - dap_result$mu) < tol & abs(log(sigma_current) - log(dap_result$sigma)) < tol){
        stop <- TRUE
      }

      t_df$lambda_loss <- c(t_df$lambda_loss, lambda_loss(dap_result, lambda_current))
      t_df$eta_loss <- c(t_df$eta_loss, eta_loss(dap_result$draws_eta, lambda_current))
      #message(sprintf("Iteration %d - dap_mu: %f ", iter_num, lambda_loss(dap_result, lambda_current), eta_loss(dap_result$draws_eta, lambda_current)))
      message(sprintf("Iteration %d - lambda loss: %f eta loss: %f normal_kl_divergence: %f", iter_num, lambda_loss(dap_result, lambda_current), eta_loss(dap_result$draws_eta, lambda_current), normal_kl_divergence(dap_result, lambda_current)))
    }

    if(stop){
      message(sprintf("Terminating self_calib on iteration %d", iter_num))
      break
    }
    mu_current <- mu_new
    sigma_current <- sigma_new
  }
  t_df <- as.data.frame(t_df)
  print(t_df)
  return(list(mu=mu_current, sigma=sigma_current, t_df=t_df))
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
