data {
  int <lower=0> J; // number of schools
  real y[J]; // estimated treatment
  real<lower=0> sigma[J]; // std of estimated effect
  real lambda_mu;
  real lambda_var;
}

transformed data {
  real lambda_sigma = sqrt(lambda_var);
}

parameters {
  vector[J] theta_trans;
  real mu; // hyper-parameter of mean
  //real log_tau;
  real<lower=0> tau;
}

transformed parameters {
  //real tau = exp(log_tau);
  vector[J] theta;
  theta = theta_trans * tau + mu;
}

model {
  theta_trans ~ normal(0, 1);
  //tau ~ cauchy(0, 5);
  tau ~ normal(lambda_mu, lambda_sigma);
  //log_tau ~ normal(lambda_mu, lambda_sigma);
  y ~ normal(theta, sigma);
  mu ~ normal(0, 5);
}
