data {
  int<lower=0> N_obs;
  int y[N_obs];

  int<lower=1> N_predictors;
  matrix[N_predictors, N_obs] x;
}

parameters {
  ordered[2] mu;

  vector[N_predictors] beta;
}

model {
  vector[N_obs] theta = inv_logit(transpose(x) * beta);


  for(n in 1:N_obs) {
    target += log_mix(theta[n],
                      poisson_log_lpmf(y[n] | mu[1]),
                      poisson_log_lpmf(y[n] | mu[2]));
  }
  target += normal_lpdf(mu | 3, 1);
  target += normal_lpdf(beta | 0, 1);
}

