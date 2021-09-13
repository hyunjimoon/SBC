data {
  // hyperparams
  int N;
  // params
  real loc_mu;
  real loc_sd;
  real scale_mu;
  real scale_sd;
  // outcome
  vector[N] y;
}

parameters {
  real loc;
  real <lower = 0> scale;
}

model {
  loc ~ normal(loc_mu, loc_sd);
  scale ~ lognormal(scale_mu, scale_sd);
  y ~ normal(loc, scale);
}
