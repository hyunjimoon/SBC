data{
  int N;
  array[N] int y;
}
parameters{
  real<lower = 0> lambda;
}
model{
  lambda ~ gamma(15, 5);
  y ~ poisson(lambda);
}
