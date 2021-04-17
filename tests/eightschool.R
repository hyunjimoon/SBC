library(cmdstanr)
library(SBC)
library(ggplot2)

J <- 8
y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)

data = list("J"=J, "y"=y, "sigma"=sigma)

n_datasets = 200
thin = 3

ncp_model = cmdstanr::cmdstan_model("eightschools_ncp.stan")
ncp_sbc_obj = SBC::SBCModel$new(name="eightschools_ncp", stan_model=ncp_model)

cp_model = cmdstanr::cmdstan_model("eightschools_cp.stan")
cp_sbc_obj = SBC::SBCModel$new(name="eightschools_cp", stan_model = cp_model)

ncp_theta_prior = ncp_sbc_obj$sample_theta_tilde_stan(list("theta", "mu", "tau"), n_datasets, data=data)
cp_theta_prior = cp_sbc_obj$sample_theta_tilde_stan(list("theta", "mu", "tau"), n_datasets, data=data)

ncp_sampled_y = ncp_sbc_obj$sample_y_tilde(ncp_theta_prior, data=data)
cp_sampled_y = cp_sbc_obj$sample_y_tilde(cp_theta_prior, data=data)

ncp_sampled_theta_post = ncp_sbc_obj$sample_theta_bar_y(ncp_sampled_y, data=data, pars=list("theta", "mu", "tau"), fit_iter=400)
cp_sampled_theta_post = cp_sbc_obj$sample_theta_bar_y(cp_sampled_y, data=data, pars=list("theta", "mu", "tau"), fit_iter=400)

ncp_rank = SBC::calculate_rank(ncp_theta_prior, ncp_sampled_theta_post, thin)
cp_rank = SBC::calculate_rank(cp_theta_prior, cp_sampled_theta_post, thin)

pars = dimnames(ncp_rank)[[2]]

for(i in 1:length(pars)){
  print(pars[[i]])
  SBC::print_summary(ncp_rank, pars[[i]], thin)
  SBC::print_summary(cp_rank, pars[[i]], thin)
  print(SBC::plot_hist(ncp_rank, pars[[i]]))
  print(SBC::plot_hist(cp_rank, pars[[i]]))
}

for(i in 1:length(pars)){
  print(pars[[i]])
  print(SBC::plot_ecdf(ncp_rank, pars[[i]]))
  print(SBC::plot_ecdf(cp_rank, pars[[i]]))
}
