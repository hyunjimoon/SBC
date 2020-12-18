test_that("calculate.rank works", {
  model_path <- file.path("sbc_array_pois.stan")
  rstan_model <- rstan::stan_model(model_path)

  rstan_sbc_obj <- SBCModel$new("poisson_rstan", rstan_model, list("lambda"=function(){as.integer(rgamma(1, shape=15, rate=5))}))
  rstan_theta_arr <- rstan_sbc_obj$sample_theta_tilde(list("lambda"), 10)

  rstan_sampled_y <- rstan_sbc_obj$sample_y_tilde(rstan_theta_arr, 25, data=list(N=25, y=as.vector(1:25)))

  rstan_sampled_y <- round(rstan_sampled_y)
  rstan_bootstrap_y <- rstan_sbc_obj$sample_bootstrap_y_tilde(rstan_sampled_y[1, ], 10)
  rstan_post <- rstan_sbc_obj$sample_theta_bar_y(rstan_bootstrap_y, pars=list("lambda"), data=list(N=25))
  expect_error(calculate.rank(rstan_theta_arr, rstan_post, 3), NA)
})
