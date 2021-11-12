data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  target += normal_lpdf(mu | 0, 1);
  target += normal_lpdf(sigma | 0, 1);
  target += normal_lpdf(y | mu, sigma);
}
