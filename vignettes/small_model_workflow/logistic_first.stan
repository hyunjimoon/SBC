data {
  int<lower=0> N_obs;
  array[N_obs] int<lower=0, upper=1> y;

  int<lower=1> N_predictors;
  matrix[N_obs, N_predictors] X;
}

parameters {
  real alpha;
  vector[N_predictors] beta;
}

model {
  for(n in 1:N_obs) {
    real mu = 0;
    for(p in 1:N_predictors) {
      mu = mu + beta[p] * X[n, p];
    }
    target += bernoulli_lpmf(y[n] | inv_logit(mu));
  }
  target += normal_lpdf(alpha | 0, 2);
  target += normal_lpdf(beta | 0, 1);
}
