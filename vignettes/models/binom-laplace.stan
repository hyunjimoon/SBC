data {
  int<lower=1> nobs;  // total number of observations
  int Y [nobs];  // outcome
  int <lower =1> nsize;
  real mu;
  real <lower = 0> sigma;
}

parameters {
  real a;
}

model {
  a ~ normal(mu, sigma);
  Y ~ binomial_logit_lpmf(nsize, a);
}
