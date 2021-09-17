data {
  int N;
  vector<lower=0>[N] y;
}

parameters {
  real<lower = 0> shape;
  real<lower = 0> scale;
}

model {
  y ~ gamma(shape, scale);
  shape ~ lognormal(0, 1);
  scale ~ lognormal(0, 1.5);
}
