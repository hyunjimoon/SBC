library(rstan)
library(cmdstanr)
library(SBC)

stan_model = "
data{
  int N;
  real Y[N];
}

parameters{
  real mu_[2];
}

model{
  Y ~ normal(mu_[1], 1);
  mu_[2] ~ normal(0, 1);
  mu_[1] ~ normal(mu_[2], 1);
}

generated quantities{
  real lambda_ = 1.0;
  vector[N] y_rep;
  for(i in 1:N){
    y_rep[i] = normal_rng(mu_[1], 1);
  }
}
"

#stan_model <- rstan::stan_model(model_code=stan_model)

cmdstan_model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stan_model))
cmdstan_sbc_obj <- SBCModel$new("poisson_cmdstan", cmdstan_model)
hypers <- list("lambda"=function(){as.integer(rgamma(1, shape=15, rate=5))})
fit <- cmdstan_model$sample(data=list(N=10, Y=c(1:10)))
cmdstan_summary <- as.data.frame(fit$summary())
row.names(cmdstan_summary) <- cmdstan_summary[, "variable"]
cmdstan_sbc_obj$infer_sequential_params(list("mu_", "lambda_"), cmdstan_summary)
#cmdstan_sbc_obj$sample_theta_tilde_stan(list("mu"), 20, list(N=10, Y=c(1,2,3,4,5,6,7,8,9,10)))
