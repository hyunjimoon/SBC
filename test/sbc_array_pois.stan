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

generated quantities{
  int y_rep[25];
  for(i in 1:25){
    y_rep[i] = poisson_rng(lambda); 
  }
  
}