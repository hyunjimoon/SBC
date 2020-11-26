library(cmdstanr)
#library(rstan)
library(posterior)
set.seed(20201111)
source(file.path(getwd(),"src/sample.R"))
model_path <- file.path(getwd(), "test/sbc_array_pois.stan")
###########################
# check cmdstanr
cmdstan_model <- cmdstan_model(model_path)

cmdstan_sbc_obj <- SBCModel$new("poisson_cmdstan", cmdstan_model, list("lambda"=function(){as.integer(rgamma(1, shape=15, rate=5))}))#, "sigma"=function(){rexp(1)}))
theta_arr <- cmdstan_sbc_obj$sample_theta_tilde(list("lambda"), 10)

cmdstan_sbc_obj$sample_y_tilde(theta_arr, 25, data=list(N=25, y=as.vector(1:25)))


###########################
# check rstan
rstan_model <- stan_model(model_path)

rstan_sbc_obj <- SBCModel$new("poisson_rstan", rstan_model, list("lambda"=function(){as.integer(rgamma(1, shape=15, rate=5))}))#, "sigma"=function(){rexp(1)}))
theta_arr <- rstan_sbc_obj$sample_theta_tilde(list("lambda"), 10)
rstan_sbc_obj$sample_y_tilde(theta_arr, 25, data=list(N=25, y=as.vector(1:25)))
