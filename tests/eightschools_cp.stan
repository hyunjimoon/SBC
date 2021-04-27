data {
  int <lower=0> J; // number of schools
  real y[J]; // estimated treatment
  real<lower=0> sigma[J]; // std of estimated effect
}
parameters {
  real theta[J]; // treatment effect in school j
  real mu; // hyper-parameter of mean
  real<lower=0> tau; // hyper-parameter of sdv
}
model {
  tau ~ cauchy(0, 5); // a non-informative prior
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
  mu ~ normal(0, 5);
}
generated quantities {
  vector[J] y_;
  vector[J] theta_;
  real mu_;
  real tau_;
  mu_ = normal_rng(0, 5);
  tau_ = fabs(cauchy_rng(0, 5));
  for(j in 1:J){
    y_[j] = normal_rng(theta[j], sigma[j]);
    theta_[j] = normal_rng(mu_, tau_);
  }
}
