data {
  int<lower=1> nobs;  // total number of observations
  int Y [nobs];  // outcome variable
  // fixed arguments for distribution
  int<lower=1> nsims; // next simulation number
  int<lower=1> nsize; // total simulation number
  // updating arguments for approximation
  vector[nsims] lambda_mu; // mean values for gmm prior
  real lambda_log_sigma; // bandwidth(sd) for gmm prior
  int<lower=1,upper=3> link;
}

parameters {
  real eta;  // population-level effects
}

transformed parameters {
  // lambda_multiple options for link functions
  real p;
  if (link == 1) {
    p = inv_logit(eta);
  } else if (link == 2) {
    p = Phi(eta);
  } else if (link == 3) {
    p = inv_cloglog(eta);
  }
}

model {
  target += binomial_logit_lpmf(Y | nsize, eta);
  // priors including constants
  target += normal_lpdf(eta | lambda_mu, exp(lambda_log_sigma)) - log(nsims);
}
