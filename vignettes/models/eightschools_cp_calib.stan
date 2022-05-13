data {
  int <lower=0> J; // number of schools
  real y[J]; // estimated treatment
  real<lower=0> sigma[J]; // std of estimated effect

  int nsims;
  vector[nsims] mm_mean;
  real mm_bandwidth;
}
parameters {
  real theta[J]; // treatment effect in school j
  real mu; // hyper-parameter of mean
  real<lower=0> tau; // hyper-parameter of sdv
}
model {
  #tau ~ cauchy(0, 5); // a non-informative prior
  target += normal_lpdf(tau | mm_mean, mm_bandwidth) - log(nsims);
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
  mu ~ normal(0, 5);
}
