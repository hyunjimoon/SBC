data {
  int<lower=1> nobs;  // total number of observations
  int Y[nobs];  // outcome
  int <lower =1> nsize;
  real mu;
  real <lower = 0> sigma;
  int<lower=1,upper=3> link;
}
parameters {
  real eta;
}
transformed parameters {
  // multiple options for link functions
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
  eta ~ normal(mu, sigma);
  Y ~ binomial(nsize, p);
}


