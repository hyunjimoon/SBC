test_that("SBCModel works for cmdstanr", {
  model_path <- file.path("sbc_array_pois.stan")

  cmdstan_model <- cmdstanr::cmdstan_model(model_path)#, force_recompile=TRUE)

  expect_error(cmdstan_sbc_obj <- SBCModel$new("poisson_cmdstan", cmdstan_model, list("lambda"=function(){as.integer(rgamma(1, shape=15, rate=5))})), NA)

  expect_error(cmd_theta_arr <- cmdstan_sbc_obj$sample_theta_tilde(list("lambda"), 10), NA)

  expect_error(cmd_sampled_y <- cmdstan_sbc_obj$sample_y_tilde(cmd_theta_arr, 25, data=list(N=25, y=as.vector(1:25))), NA)
  expect_error(cmd_bootstrap_y <- cmdstan_sbc_obj$sample_bootstrap_y_tilde(cmd_sampled_y[1, ], 10), NA)
  expect_error(cmdstan_sbc_obj$sample_theta_bar_y(cmd_bootstrap_y, pars=list("lambda", "lp__"), data=list(N=25)), NA)
})

test_that("SBCModel works for rstan", {
  model_path <- file.path("sbc_array_pois.stan")
  rstan_model <- rstan::stan_model(model_path)

  expect_error(rstan_sbc_obj <- SBCModel$new("poisson_rstan", rstan_model, list("lambda"=function(){as.integer(rgamma(1, shape=15, rate=5))})), NA)
  expect_error(rstan_theta_arr <- rstan_sbc_obj$sample_theta_tilde(list("lambda"), 10), NA)

  expect_error(rstan_sampled_y <- rstan_sbc_obj$sample_y_tilde(rstan_theta_arr, 25, data=list(N=25, y=as.vector(1:25))), NA)

  rstan_sampled_y <- round(rstan_sampled_y)
  expect_error(rstan_bootstrap_y <- rstan_sbc_obj$sample_bootstrap_y_tilde(rstan_sampled_y[1, ], 10), NA)
  expect_error(rstan_post <- rstan_sbc_obj$sample_theta_bar_y(rstan_bootstrap_y, pars=list("lambda"), data=list(N=25)), NA)

})
