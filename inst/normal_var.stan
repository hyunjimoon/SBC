data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  y ~ normal(mu, sigma^2);
  mu ~ normal(0, 1);
  sigma ~ normal(0, 1);
}

