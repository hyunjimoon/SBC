data {
  int<lower=1> nobs;  // total number of observations
  int Y [nobs];  // outcome variable
  // fixed arguments for distribution
  int<lower=1> nsims; // next simulation number
  int<lower=1> nsize; // total simulation number
  // updating arguments for approximation
  vector[nsims] mixture_means; // mean values for gmm prior
  real mixture_sds; // bandwidth(sd) for gmm prior
}

parameters {
  real a;  // population-level effects
}

model {
  target += binomial_logit_lpmf(Y | nsize, a);
  // priors including constants
  target += normal_lpdf(a | mixture_means, mixture_sds) - log(nsims);
}
