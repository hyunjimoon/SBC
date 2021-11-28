data {
  int <lower=0> J; // number of schools
  real y[J]; // estimated treatment
  real<lower=0> sigma[J]; // std of estimated effect
  real lambda_mu;
  real lambda_log_sigma;
}
parameters {
  real theta[J]; // treatment effect in school j
  real mu; // hyper-parameter of mean
  real log_tau;
  //real<lower=0> tau; // hyper-parameter of sdv
}
//transformed parameters {
//  real tau = exp(-0.5);
//}

model {
  //tau ~ cauchy(0, 5);
  log_tau ~ normal(lambda_mu, exp(lambda_log_sigma));
  theta ~ normal(mu, exp(log_tau));
  y ~ normal(theta, sigma);
  //mu ~ normal(lambda_mu, exp(lambda_log_sigma));
  mu ~ normal(0, 5);
}
