library(cmdstanr)
library(SBC)
library(posterior)

matrix_test_code <- "
data {
  int J; // estimated treatment
}
parameters {
}

transformed parameters {
  matrix[2, 3] A;
  vector[5] B;
  for(i in 1:2){
    for(j in 1:3){
      A[i, j] = 10 * i + j;
    }
  }
  for(i in 1:5){
    B[i] = i;
  }
}

model {
}
generated quantities {
  matrix[2, 3] A_;
  vector[5] B_;
  for(i in 1:2){
    for(j in 1:3){
      A_[i, j] = 10 * i + j;
    }
  }
  for(i in 1:5){
    B_[i] = i;
  }
}

"

model <- cmdstan_model(write_stan_file(matrix_test_code, basename = "SBC_matrix_test"))
J <- 8

data = list("J"=J)

n_datasets = 2
thin = 3

sbc_obj = SBC::SBCModel$new(name="test", stan_model=model)

res <- model$sample(data=data, fixed_param = TRUE)
sample_summary <- as.data.frame(res$summary())
row.names(sample_summary) <- sample_summary[, "variable"]
sbc_obj$infer_sequential_params(list("A", "B"), sample_summary, return_dim_info = TRUE)


theta_prior = sbc_obj$sample_theta_tilde_stan(list("A", "B"), n_datasets, data=data)

sampled_y = sbc_obj$sample_y_tilde(theta_prior, data=data)

#ncp_sampled_theta_post = ncp_sbc_obj$sample_theta_bar_y(ncp_sampled_y, data=data, pars=list("theta", "mu", "tau"), fit_iter=400)

#ncp_rank = SBC::calculate_rank(ncp_theta_prior, ncp_sampled_theta_post, thin)

# pars = dimnames(ncp_rank)[[2]]
#
# for(i in 1:length(pars)){
#   print(pars[[i]])
#   SBC::print_summary(ncp_rank, pars[[i]], thin)
#   SBC::print_summary(cp_rank, pars[[i]], thin)
#   print(SBC::plot_hist(ncp_rank, pars[[i]]))
#   print(SBC::plot_hist(cp_rank, pars[[i]]))
# }
#
# for(i in 1:length(pars)){
#   print(pars[[i]])
#   print(SBC::plot_ecdf(ncp_rank, pars[[i]]))
#   print(SBC::plot_ecdf(cp_rank, pars[[i]]))
# }

