data {
  int<lower=1> nobs;  // total number of observations
  vector[nobs] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[nobs, K] X;  // population-level design matrix
  //real<lower=0> shape;  // shape parameter
  int dist_types[2];
  real lambda_arg1[2];
  real lambda_arg2[2];



}
parameters {
  real a;  // population-level effects
  vector[K] b;
  real<lower=1e-100> shape;
}

model {
  // initialize linear predictor term
  vector[nobs] logmu = a + X * b;
  vector[nobs] mu;
  for (n in 1:nobs) {
    // apply the inverse link function
    mu[n] = exp(logmu[n]);
  }

  for(n in 1:nobs){
    target += gamma_lpdf(Y[n] | shape, shape / mu[n]);
  }

  // priors including constants
  //target += gamma_lpdf(shape | lambda_alpha, lambda_beta);
  //target += normal_lpdf(a | 2, 5);
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


  if (dist_types[1] == 1){ // shape
    target += normal_lpdf(shape | lambda_arg1[1], lambda_arg2[1]);
  }
  else if(dist_types[2] == 2){
    target += gamma_lpdf(shape | lambda_arg1[1], lambda_arg2[1]);
  }


  if (dist_types[2] == 1){ // a
    target += normal_lpdf(a | lambda_arg1[2], lambda_arg2[2]);
  }
  else if(dist_types[2] == 2){
    target += normal_lpdf(a | lambda_arg1[2], lambda_arg2[2]);
  }
}
