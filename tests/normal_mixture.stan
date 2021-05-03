// normal mixture, unknown proportion and means, known variance
// p(y|mu,theta) = theta * Normal(y|mu[1],1) + (1-theta) * Normal(y|mu[2],1);

data {
  int<lower=0>  N;
  real y[N];
}
parameters {
  real<lower=0,upper=1> theta;
  real mu[2];
}
model {
  theta ~ uniform(0,1); // equivalently, ~ beta(1,1);
  for (k in 1:2)
    mu[k] ~ normal(0,10);
  for (n in 1:N)
    target += log_mix(theta, normal_lpdf(y[n]|mu[1],1.0), normal_lpdf(y[n]|mu[2],1.0));
}

generated quantities {
  real theta_;
  real mu_[2];
  vector[N] y_;
  theta_ = uniform_rng(0, 1);
  for (k in 1:2){
    mu_[k] = normal_rng(0,10);
  }
  for(j in 1:N){
    y_[j] = normal_rng(mu[1], 1.0) * theta + normal_rng(mu[2], 1.0) * (1 - theta);
  }


}
