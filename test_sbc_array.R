library(cmdstanr)
library(posterior)
set.seed(20201111)
source("sbc_array.R")
n_iter <- 100
n_pars <- 1
n_sample <- 100

prior_lambdas <- abs(rnorm(n_iter, 0, 2))#seq(1, n_iter)
prior_arr <- array(dim=c(n_iter, n_pars))
dimnames(prior_arr)[2] <- list("lambda")
prior_arr[, "lambda"] <- prior_lambdas
model_path <- "sbc_array_pois.stan"
model <- cmdstan_model(model_path)

draw_list <- array(dim=c(n_iter, n_pars, n_sample))
dimnames(draw_list)[2] <- list("lambda")
for(i in 1:n_iter){
  prior_samples <- rpois(n_sample, prior_lambdas[i])
  model_fit <- model$sample(data=list(N=n_sample, y=prior_samples), iter_warmup = n_sample, iter_sampling = n_sample,chains=1, parallel_chains=1, save_warmup = FALSE)#, output_dir=file.path(getwd(), "/output_dir"))
  res <- as_draws_array(model_fit$draws(variables=c("lambda")))
  draw_list[i, , ] <- res
}
print(prior_arr)
sbc <- sbc.rank(prior_arr, draw_list, list("lambda"), 3)
sbc.plot(sbc, "lambda", 3)
