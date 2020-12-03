library(cmdstanr)
library(rstan)
library(posterior)
set.seed(20201111)
source(file.path(getwd(),"src/sample.R"))
model_path <- file.path(getwd(), "test/sbc_array_pois.stan")
###########################
# check cmdstanr
cmdstan_model <- cmdstan_model(model_path)#, force_recompile=TRUE)

cmdstan_sbc_obj <- SBCModel$new("poisson_cmdstan", cmdstan_model, list("lambda"=function(){as.integer(rgamma(1, shape=15, rate=5))}))#, "sigma"=function(){rexp(1)}))
cmd_theta_arr <- cmdstan_sbc_obj$sample_theta_tilde(list("lambda"), 10)

cmd_sampled_y <- cmdstan_sbc_obj$sample_y_tilde(cmd_theta_arr, 25, data=list(N=25, y=as.vector(1:25)))
cmd_bootstrap_y <- cmdstan_sbc_obj$sample_bootstrap_y_tilde(cmd_sampled_y[1, ], 10)
cmd_post <- cmdstan_sbc_obj$sample_theta_bar_y(cmd_bootstrap_y, pars=list("lambda", "lp__"), data=list(N=25))
print(dim(cmd_post))
print(cmd_theta_arr[1, "lambda"])
print(mean(cmd_post[, "lambda", ]))
###########################
# check rstan
rstan_model <- stan_model(model_path)

rstan_sbc_obj <- SBCModel$new("poisson_rstan", rstan_model, list("lambda"=function(){as.integer(rgamma(1, shape=15, rate=5))}))#, "sigma"=function(){rexp(1)}))
rstan_theta_arr <- rstan_sbc_obj$sample_theta_tilde(list("lambda"), 10)
rstan_sampled_y <- rstan_sbc_obj$sample_y_tilde(rstan_theta_arr, 25, data=list(N=25, y=as.vector(1:25)))
rstan_sampled_y <- round(rstan_sampled_y)  # rstan returns double values(float components are ignoreable, < 1e-4)
rstan_bootstrap_y <- rstan_sbc_obj$sample_bootstrap_y_tilde(rstan_sampled_y[1, ], 10)
rstan_post <- rstan_sbc_obj$sample_theta_bar_y(rstan_bootstrap_y, pars=list("lambda"), data=list(N=25))
print(dim(rstan_post))
print(rstan_theta_arr[1, "lambda"])
print(mean(rstan_post[, "lambda", ]))

################################
source(file.path(getwd(),"src/sbc_array.R"))
#sbc.obj <- new("SBCData", prior=cmd_theta_arr, posterior=cmd_post[, "lambda", ], model.name="poisson")
rstan_theta_arr[2:10, ] <- rstan_theta_arr[1, ]
sbc.obj <- new("SBCData", prior=rstan_theta_arr, posterior=rstan_post, model.name="poisson")
rank <- sbc.rank(sbc.obj, 3)
sbc.plot.hist(rank, "lambda", 3)
sbc.plot.ecdf(rank, "lambda")
sbc.plot.ecdf.diff(rank, "lambda")  # this plot not looking good :(
