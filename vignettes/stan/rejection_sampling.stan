data {
   int<lower=0> N;
   real y[N];
}

parameters {
   real mu;
}

model {
   mu ~ normal(0, 2);
   y ~ normal(mu, 1);
}
