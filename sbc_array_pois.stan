data{
  int N; 
  int y[N];
}
parameters{
  real<lower = 0> lambda;
}
model{
  y ~ poisson(lambda);
}
