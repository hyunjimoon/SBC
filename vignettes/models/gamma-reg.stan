data {
  int<lower=1> nobs;  // total number of observations
  vector[nobs] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[nobs, K] X;  // population-level design matrix
  real<lower=0> shape;  // shape parameter
  real lambda_mu;
  real lambda_log_sigma;
}
parameters {
  real eta;  // population-level effects
  vector[K] b;
}

model {
  // initialize linear predictor term
  vector[nobs] logmu = eta + X * b;
  vector[nobs] mu;
  for (n in 1:nobs) {
    // apply the inverse link function
    mu[n] = shape * exp(-(logmu[n]));
  }
  target += gamma_lpdf(Y | shape, mu);

  // priors including constants
  target += normal_lpdf(eta| lambda_mu, exp(lambda_log_sigma));
  target += normal_lpdf(b[1] | 0, 1);
  target += normal_lpdf(b[2] | 0, 1);
  target += normal_lpdf(b[3] | 0, 1);
  target += normal_lpdf(b[4] | 0, 1);
  target += normal_lpdf(b[5] | 0, 1);
  target += normal_lpdf(b[6] | 0, 1);
  target += normal_lpdf(b[7] | 0, 1);
  target += normal_lpdf(b[8] | 0, 1);
  target += normal_lpdf(b[9] | 0, 1);
  target += normal_lpdf(b[10] | 0, 1);
  target += normal_lpdf(b[11] | 0, 1);
  target += normal_lpdf(b[12] | 0, 1);
  target += normal_lpdf(b[13] | 0, 1);
  target += normal_lpdf(b[14] | 0, 1);
  target += normal_lpdf(b[15] | 0, 1);
}
