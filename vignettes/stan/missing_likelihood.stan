data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real mu;
}

model {
  target += normal_lpdf(mu | 0, 1);
}
