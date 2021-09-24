
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  real<lower=0> shape;  // shape parameter

  int<lower=1> n_calib_parameters;
  int SM;
  matrix[n_calib_parameters, SM] mm_mean; // mean values for gmm prior
  vector[n_calib_parameters] mm_bandwidth; // bandwidth(sd) for gmm prior
}
parameters {
  real a;  // population-level effects
  vector[K] b;
}

model {
  // initialize linear predictor term
  vector[N] mu = a + X * b;
  for (n in 1:N) {
    // apply the inverse link function
    mu[n] = shape * exp(-(mu[n]));
  }
  target += gamma_lpdf(Y | shape, mu);

  // priors including constants
  //target += normal_lpdf(a| a_loc, a_scale);

  target += normal_lpdf(a | row(mm_mean, 1), mm_bandwidth[1]) - log(SM);

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
