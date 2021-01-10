library(SBC)
library(cmdstanr)
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

# set hyperprior and iteration numbers
n_datasets = 10 # total number of simulated parameter sets from theta
thin = 3

data <- list(n_datasets=n_datasets, y=as.vector(1:n_datasets))

model <- cmdstan_model(write_stan_file(stan_code))
# alternatively, use rstan:
# model <- stan_model(model_code = stan_code)

sbc_obj <- SBCModel$new(name = "poisson", stan_model = model)


hyperpriors <- list("lambda"=function(){rgamma(1, shape=15, scale=0.2)})
theta_prior <- sbc_obj$sample_theta_tilde_stan(list("lambda"), n_datasets, data=data)
# or alternatively, use the hyperpriors list:
# theta_prior <- sbc_obj$sample_theta_tilde_stan(list("lambda), n_datasets, hyperpriors)

sampled_y <- sbc_obj$sample_y_tilde(theta_prior , data=data)

theta_post <- sbc_obj$sample_theta_bar_y(sampled_y, data=data, pars=list("lambda"), fit_iter = 200)

rank <- calculate_rank(theta_prior, theta_post, thin = thin)
bootstrap_y <- sbc_obj$sample_bootstrap_y_tilde(sampled_y[1, ], 10)  # create 10 sample vectors from a single sample vector
