library(cmdstanr)
library(rstan)
library(posterior)
set.seed(20201111)
source(file.path(getwd(),"src/sample.R"))
model_path <- file.path(getwd(), "test/sbc_array_pois.stan")
###########################
# check cmdstanr
cmdstan_model <- cmdstan_model(model_path)

cmdstan_sbc_obj <- SBCModel$new("poisson_cmdstan", cmdstan_model, list("lambda"=function(){as.integer(rgamma(1, shape=15, rate=5))}))#, "sigma"=function(){rexp(1)}))
theta_arr <- cmdstan_sbc_obj$sample_theta_tilde(list("lambda"), 1)

sampled_y <- cmdstan_sbc_obj$sample_y_tilde(theta_arr, 25, data=list(N=25, y=as.vector(1:25)))

post <- cmdstan_sbc_obj$approx_theta_bar_y(sampled_y, pars=list("lambda", "lp__"), data=list(N=25))
print(dim(post))
print(theta_arr[1, "lambda"])
print(mean(post[, "lambda", ]))
###########################
# check rstan
rstan_model <- stan_model(model_path)

rstan_sbc_obj <- SBCModel$new("poisson_rstan", rstan_model, list("lambda"=function(){as.integer(rgamma(1, shape=15, rate=5))}))#, "sigma"=function(){rexp(1)}))
theta_arr <- rstan_sbc_obj$sample_theta_tilde(list("lambda"), 100)
sampled_y <- rstan_sbc_obj$sample_y_tilde(theta_arr, 25, data=list(N=25, y=as.vector(1:25)))
sampled_y <- round(sampled_y)  # rstan returns double values(float components are ignoreable, < 1e-4)
post_samples <- rstan_sbc_obj$approx_theta_bar_y(sampled_y, pars=list("lambda"), data=list(N=25))
print(dim(post_samples))
print(theta_arr[1, "lambda"])
print(mean(post_samples[, "lambda", ]))

################################
source(file.path(getwd(),"src/sbc_array.R"))
sbc.obj <- new("SBCData", prior=theta_arr, posterior=post_samples, model.name="poisson")
rank <- sbc.rank(sbc.obj, 3)
sbc.plot.hist(rank, "lambda", 3)
sbc.plot.ecdf(rank, "lambda")
sbc.plot.ecdf.diff(rank, "lambda")  # this plot not looking good :(
