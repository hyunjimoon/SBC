data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of predictors
  matrix[N, K] x;   // predictor matrix
  vector[N] y;      // outcome vector
}
parameters {
  real alpha;           // intercept
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  // error scale
}
model {
  vector[N] mu;
  for(i in 1:N) {
    mu[i] = alpha;
    for(j in 1:K) {
      mu[i] += beta[j] * x[j, j];
    }
  }
  y ~ normal(mu, sigma);  // likelihood
  alpha ~ normal(0, 5);
  beta ~ normal(0, 1);
  sigma ~ normal(0, 2);
}
