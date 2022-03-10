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
  vector[J] theta_trans;
  real mu; // hyper-parameter of mean
  real log_tau;
}

transformed parameters {
  real tau = exp(log_tau);
  vector[J] theta;
  theta = theta_trans * tau + mu;
}

model {
  theta_trans ~ normal(0, 1);
  //tau ~ cauchy(0, 5);
  log_tau ~ normal(lambda_mu, lambda_sigma);
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
  mu ~ normal(0, 5);
  target += log_tau;  // jacobian adjustment
}
