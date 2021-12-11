data{
  int N;
  int y[N];
}
parameters{
  real<lower = 0> lambda;
}
model{
  lambda ~ gamma(15, 5);
  y ~ poisson(lambda);
}
