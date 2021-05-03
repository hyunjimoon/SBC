library(cmdstanr)
library(SBC)
library(ggplot2)


data <- rjson::fromJSON(file="tests/normal_2.json")

n_datasets = 200
thin = 3

model = cmdstanr::cmdstan_model("tests/normal_mixture.stan")
sbc_obj = SBC::SBCModel$new(name="normal_mix_2", stan_model=model)

theta_prior = sbc_obj$sample_theta_tilde_stan(list("theta", "mu"), n_datasets, data=data)

sampled_y = sbc_obj$sample_y_tilde(theta_prior, data=data)

sampled_theta_post = sbc_obj$sample_theta_bar_y(sampled_y, data=data, pars=list("theta", "mu"), fit_iter=25*3)

rank = SBC::calculate_rank(theta_prior, sampled_theta_post, thin)

pars = dimnames(rank)[[2]]

for(i in 1:length(pars)){
  print(pars[[i]])
  SBC::print_summary(rank, pars[[i]], thin)
  print(SBC::plot_hist(rank, pars[[i]]))
}

for(i in 1:length(pars)){
  print(pars[[i]])
  print(SBC::plot_ecdf(rank, pars[[i]]))
}
