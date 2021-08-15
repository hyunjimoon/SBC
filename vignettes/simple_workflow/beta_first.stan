data {
  int<lower=0> N_obs;
  vector<lower=0, upper=1>[N_obs] y;

  int<lower=1> N_predictors;
  matrix[N_predictors, N_obs] x;
}

parameters {
  matrix[2, N_predictors] beta;
}

model {
  matrix[2, N_obs] linpred = beta * x;
  target += beta_lpdf(y | exp(linpred[1,]), exp(linpred[2,]));
  target += normal_lpdf(to_vector(beta) | 0, 1);
}
