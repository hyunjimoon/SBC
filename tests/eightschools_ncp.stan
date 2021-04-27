data {
  int <lower=0> J; // number of schools
  real y[J]; // estimated treatment
  real<lower=0> sigma[J]; // std of estimated effect
}
parameters {
  vector[J] theta_trans; // transformation of theta
  real mu; // hyper-parameter of mean
  real<lower=0> tau; // hyper-parameter of sd
}
transformed parameters{
  vector[J] theta;
  // original theta
  theta=theta_trans*tau+mu;
}
model {
  theta_trans ~ normal (0,1);
  y ~ normal(theta , sigma);
  mu ~ normal(0, 5); // a non-informative prior
  tau ~ cauchy(0, 5);
}

generated quantities {
  vector[J] y_;
  vector[J] theta_trans_;
  vector[J] theta_;
  real mu_;
  real tau_;
  mu_ = normal_rng(0, 5);
  tau_ = fabs(cauchy_rng(0, 5));
  for(j in 1:J){
    y_[j] = normal_rng(theta[j], sigma[j]);
    theta_trans_[j] = normal_rng(0, 1);
    theta_[j] = theta_trans_[j] * tau_ + mu_;
  }
}
