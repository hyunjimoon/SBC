data {
  int N; // Number of observations
  array[N] int y;
}
parameters {
  // Parameters of measurement model
  ordered[2] log_mu;

  // Initial state
  simplex[2] rho;

  // Rows of the transition matrix
  simplex[2] t1;
  simplex[2] t2;
}

model {

  matrix[2, 2] Gamma;
  matrix[2, N] log_omega;

  // Build the transition matrix
  Gamma[1, : ] = t1';
  Gamma[2, : ] = t2';

  // Compute the log likelihoods in each possible state
  for (n in 1 : N) {
    // The observation model could change with n, or vary in a number of
    //  different ways (which is why log_omega is passed in as an argument)
    log_omega[1, n] = poisson_log_lpmf(y[n] | log_mu[1]);
    log_omega[2, n] = poisson_log_lpmf(y[n] | log_mu[2]);
  }

  log_mu[1] ~ normal(-2, 1);
  log_mu[2] ~ normal(2, 1);

  // Initial state - we're quite sure we started with the source working
  rho ~ dirichlet([1, 10]);

  t1 ~ dirichlet([3, 3]);
  t2 ~ dirichlet([3, 3]);

  target += hmm_marginal(log_omega, Gamma, rho);
}


generated quantities {
  positive_ordered[2] mu = exp(log_mu);
}
