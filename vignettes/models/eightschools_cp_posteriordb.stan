data {
  int <lower=0> J; // number of schools
  real y[J]; // estimated treatment
  real<lower=0> sigma[J]; // std of estimated effect
  real lambda_mu;
  real lambda_log_var;
}

transformed data {
  real lambda_sigma = sqrt(exp(lambda_log_var));
}

parameters {
  real theta[J]; // treatment effect in school j
  real mu; // hyper-parameter of mean
  real log_tau;
}

transformed parameters {
  real tau = exp(log_tau);
}

model {
  //tau ~ cauchy(0, 5);
  log_tau ~ normal(lambda_mu, lambda_sigma);
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
  mu ~ normal(0, 5);
  target += log_tau;  // jacobian adjustment
}
