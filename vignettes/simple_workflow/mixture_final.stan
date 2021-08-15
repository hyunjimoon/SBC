data {
  int<lower=0> N;
  int<lower=2> N_components;
  int y[N];
}

parameters {
  ordered[N_components] mu;
  simplex[N_components] theta;
}

model {
  for(n in 1:N) {
    vector[N_components] component_lp;
    for(c in 1:N_components) {
      component_lp[c] = poisson_log_lpmf(y[n] | mu[c]);
    }
    target += log_sum_exp(theta .* component_lp);
  }
  target += normal_lpdf(mu | 5, 2.5);
  target += dirichlet_lpdf(theta | rep_vector(3, N_components));
}

