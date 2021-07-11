## ----setup, include = FALSE---------------------------------------------------
library(SBC); library(cmdstanr)

## ---- include = TRUE----------------------------------------------------------
stan_code <- "
data{
  int n_datasets;
  int y[n_datasets];
}
parameters{
  real<lower = 0> lambda;
}
model{
  lambda ~ gamma(15, 5);
  y ~ poisson(lambda);
}

generated quantities{
  real lambda_ = gamma_rng(15, 5);  // optional
  int y_[n_datasets];
  real log_lik[n_datasets];
  for (i in 1:n_datasets){
    y_[i] = poisson_rng(lambda); 
    log_lik[i] = poisson_lpmf(y_[i] | lambda);
  }
}
"

## -----------------------------------------------------------------------------
# set hyperprior and iteration numbers
n_datasets = 100 # total number of simulated parameter sets from theta
thin = 3

data <- list(n_datasets=n_datasets, y=as.vector(1:n_datasets))

model <- cmdstan_model(write_stan_file(stan_code))
# alternatively, use rstan:
# model <- stan_model(model_code = stan_code)

sbc_obj <- SBCModel$new(name = "poisson", stan_model = model)

## -----------------------------------------------------------------------------
hyperpriors <- list("lambda"=function(){rgamma(1, shape=15, scale=0.2)})

## ---- results=FALSE-----------------------------------------------------------
theta_prior <- sbc_obj$sample_theta_tilde_stan(list("lambda"), n_datasets, data = data)

# or alternatively, use the hyperpriors list:
# theta_prior <- sbc_obj$sample_theta_tilde(list("lambda), n_datasets, hyperpriors)

## ---- results=FALSE-----------------------------------------------------------
sampled_y <- sbc_obj$sample_y_tilde(theta_prior , data=data)

## ---- results=FALSE-----------------------------------------------------------
# retrieve parameters by refitting on simulated y data
theta_post <- sbc_obj$sample_theta_bar_y(sampled_y, data=data, pars=list("lambda"), fit_iter = 200)

## -----------------------------------------------------------------------------
rank <- calculate_rank(theta_prior, theta_post, thin = thin)  # thinning factor of 3
plot_hist(rank, "lambda")
plot_ecdf(rank, "lambda")

## ---- results=FALSE-----------------------------------------------------------
model <- cmdstan_model(write_stan_file(stan_code))
# alternatively, use rstan:
# model <- stan_model(model_code = stan_code)

single_sbc_obj <- SBCModel$new(name = "poisson", stan_model = model)
easy_ranks <- single_sbc_obj$sample_all(priors=list(), pars=list("lambda"), n_iters=20, n_fits=200, data=data, thin=3)

## -----------------------------------------------------------------------------
plot_hist(easy_ranks, "lambda")
plot_ecdf(easy_ranks, "lambda")

## ---- results=FALSE-----------------------------------------------------------
bootstrap_y <- sbc_obj$sample_bootstrap_y_tilde(sampled_y[1, ], 10)  # create 10 sample vectors from a single sample vector
post <- sbc_obj$sample_theta_bar_y(bootstrap_y, data=data, pars=list("lambda"), fit_iter=200)
# Do the rank magic here

