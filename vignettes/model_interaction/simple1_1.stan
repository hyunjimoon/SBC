data {
  int N;
  vector[N] y;
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  sigma ~ lognormal(0, 1);
  mu ~ normal(0, 1);
  y ~ normal(mu, sigma);
}
