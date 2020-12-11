library(cmdstanr)
library(posterior)
set.seed(20201111)
source(file.path(getwd(),"src/sbc_array.R"))
n_iter <- 300
n_pars <- 1
n_sample <- 200

prior_lambdas <- abs(rnorm(n_iter, 0, 2))#seq(1, n_iter)
prior_arr <- array(dim=c(n_iter, n_pars))
dimnames(prior_arr)[2] <- list("lambda")
prior_arr[, "lambda"] <- prior_lambdas
model_path <- file.path(getwd(), "test/sbc_array_pois.stan")
model <- cmdstan_model(model_path)
print(paste(typeof(model), class(model)))
draw_list <- array(dim=c(n_iter, n_pars, n_sample))
dimnames(draw_list)[2] <- list("lambda")
for(i in 1:n_iter){
  prior_samples <- rpois(n_sample, prior_lambdas[i])
  model_fit <- model$sample(data=list(N=n_sample, y=prior_samples), iter_warmup = n_sample, iter_sampling = n_sample,chains=1, parallel_chains=1, save_warmup = FALSE, thin=1, refresh=0)#, output_dir=file.path(getwd(), "/output_dir"))
  res <- as_draws_array(model_fit$draws(variables=c("lambda")))
  draw_list[i, , ] <- res
}

sbc.obj <- new("SBCData", prior=prior_arr, posterior=draw_list, model.name="poisson")

rank <- sbc.rank(sbc.obj, 3)
sbc.summary(rank, "lambda")
sbc.plot.hist(rank, "lambda", 3)
sbc.plot.ecdf(rank, "lambda")
sbc.plot.ecdf.diff(rank, "lambda")
