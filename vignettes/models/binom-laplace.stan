data {
  int<lower=1> nobs;  // total number of observations
  int Y[nobs];  // outcome
  int <lower =1> nsize;
  real lambda_mu;
  real lambda_log_sigma;
  int<lower=1,upper=3> link;
  int eta_dist_type;  // 1 == normal, 2 == gamma

}
parameters {
  real eta;
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
  eta ~ normal(lambda_mu, exp(lambda_log_sigma));
  Y ~ binomial(nsize, p);
}
