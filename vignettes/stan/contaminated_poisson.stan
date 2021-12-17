data {
  int N; // Number of observations
  array[N] int y;
}
parameters {
  // Parameters of measurement model
  real<lower=0> mu_signal;
  real<lower=0> mu_background;

  real<lower=0, upper=1> p_signal;
}

model {

  for (n in 1 : N) {
    target += log_mix(p_signal ,
      poisson_lpmf(y[n] | mu_background + mu_signal),
      poisson_lpmf(y[n] | mu_background));
  }

  mu_background ~ lognormal(-2, 0.2);
  mu_signal ~ lognormal(2, 1);

  p_signal ~ beta(2, 2);
}
