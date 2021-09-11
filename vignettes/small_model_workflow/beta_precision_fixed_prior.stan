data {
  int<lower=0> N_obs;
  vector<lower=0, upper=1>[N_obs] y;

  int<lower=1> N_predictors;
  matrix[N_predictors, N_obs] x;
}

parameters {
  vector[N_predictors] beta;
  real<lower=0> phi;
}

model {
  vector[N_obs] mu = inv_logit(transpose(x) * beta);

  target += beta_lpdf(y | mu * phi, (1 - mu) * phi);
  target += normal_lpdf(beta | 0, 1);
  target += lognormal_lpdf(phi | 3, 1);
}
