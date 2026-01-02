data {
  int<lower = 1> N;
  array[N] int<lower = 0, upper = 1> y;
  vector[N] x;
  int<lower = 1> C;
  array[N] int<lower = 1, upper = C> clutch;
}

parameters {
  real alpha0_raw;
  real alpha1_raw;
  vector[C - 1] b_raw;
}

transformed parameters {
  vector[C] b = append_row(b_raw, 0);
  real alpha0 = sqrt(10) * alpha0_raw;
  real alpha1 = sqrt(10) * alpha1_raw;
}

model {
  // priors

  target += normal_lpdf(alpha0_raw | 0, 1);
  target += normal_lpdf(alpha1_raw | 0, 1);

  // fixed effects
  target += normal_lpdf(b_raw | 0, 1);

  // likelihood
  for (i in 1:N) {
    target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1 * x[i] + b[clutch[i]]));
  }
}
