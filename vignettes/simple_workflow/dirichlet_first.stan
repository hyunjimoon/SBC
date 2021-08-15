data {
  int<lower=0> N_obs;
  int<lower=2> N_components;
  matrix<lower=0, upper=1>[N_components, N_obs] y;

  int<lower=1> N_predictors;
  matrix[N_predictors, N_obs] x;
}

parameters {
  matrix[N_components, N_predictors] beta;
  real<lower=0> precision;
}

model {
  matrix[N_components, N_obs] linpred = beta * x;
  for(n in 1:N_obs) {
    target += dirichlet_lpdf(y[, n] | precision * softmax(linpred[,n]));
  }

  target += lognormal_lpdf(precision | 2, 1);
  target += normal_lpdf(to_vector(beta) | 0, 2);
}
