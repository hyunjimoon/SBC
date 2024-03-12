data {
   int<lower=0> N;
   array[N] real y;
}

parameters {
   real mu;
}

model {
   mu ~ normal(0, 2);
   y ~ normal(mu, 1);
}
