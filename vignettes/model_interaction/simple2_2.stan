data {
  int N;
  vector[N] y;
}
parameters {
}
model {
  y ~ normal(0, 1);
}
