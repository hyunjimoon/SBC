data {
  int<lower=0> N;
  vector[N] y;
}

transformed data {
  int N2 = N / 2 + 1;
}

parameters {
  real mu;
}

model {
  target += normal_lpdf(mu | 0, 1);
  for(n in 1:N2) {
    target += normal_lpdf(y[n] | mu, 1);
  }
}
