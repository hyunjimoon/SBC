data {
  int N;
  vector[N] y;
}
parameters {
  real<lower=0> sigma;
}
transformed parameters{
  real<lower=0> mu;
  mu = 0;
}
model {
  sigma ~ lognormal(0, 1);
  y ~ normal(0, sigma);
}
