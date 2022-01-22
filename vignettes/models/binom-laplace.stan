data {
  int<lower=1> nobs;  // total number of observations
  int Y[nobs];  // outcome
  int <lower =1> nsize;
  int<lower=1,upper=3> link;
  int dist_types;
  real lambda_arg1;
  real lambda_arg2;

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
  if(dist_types == 1){
    #eta ~ normal(lambda_arg1, lambda_arg2);
    eta ~ gamma(lambda_arg1, lambda_arg2);
  }
  else if(dist_types == 2){
    eta ~ gamma(lambda_arg1, lambda_arg2);
  }
  else{
    reject("non-supportive distribution type");
  }


  Y ~ binomial(nsize, p);
}
