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
  vector[C] b_raw;
  //real<lower = 0> sigma2;
  real<lower = 0> sigma;
}
transformed parameters {
  vector[C] b;
  //real<lower = 0> sigma = sqrt(sigma2);
  real alpha0 = sqrt(10.0) * alpha0_raw;
  real alpha1 = sqrt(10.0) * alpha1_raw;
  b = sigma * b_raw;
}
model {
  // priors
  //target += - 2 * log1p(sigma2) -log(pi()/4); // p(sigma2) = 1 / (1 + sigma2) ^ 2
  target += normal_lpdf(sigma | 0, 1) + log(2);
  target += normal_lpdf(alpha0_raw | 0, 1);
  target += normal_lpdf(alpha1_raw | 0, 1);

  // random effects
  target += normal_lpdf(b_raw | 0, 1);

  // likelihood
  for (i in 1:N) {
    target += bernoulli_logit_lpmf(y[i] | alpha0 + alpha1 * x[i] + b[clutch[i]]);
    //target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1 * x[i] + b[clutch[i]]));
  }
}
