data {
  int N;
  vector[N] y;
}
parameters {
  real mu;
  real <lower =0> sigma;
}
model {
  mu ~ normal(0, 1);
  sigma ~ exponential(1);
  y ~ normal(mu, sigma);
}
