library(cmdstanr)

# prior sample: N samples, for each iteration
# posterior sample: N * Y samples, Y sample in each iteration
# rank = sum(prior < posterior)

# sbc.cmdstanr_posterior2arr(cmdstanfit, pars){
#   draws <- cmdstanfit$draws(variables=pars)
# }

sbc.rank <- function(prior_array, posterior_array, pars, thin){
  # n_iter: number of iterations
  # n_sample: number of posterior samples per iteration
  # n_pars: number of parameters of interest
  # prior array dimensions: (n_iter, n_pars)
  # posterior array dimensions: (n_iter, n_pars, n_sample)
  # return array dimensions: (n_iter, n_pars)
  prior_dims = dim(prior_array)
  posterior_dims = dim(posterior_array)
  if (prior_dims[1] != posterior_dims[1] || prior_dims[2] != posterior_dims[2]){
    stop("Dimension mismatch error!\n dim(prior_array)[1] == dim(posterior_array)[1] && dim(prior_array)[2] == dim(posterior_array)[2] must be satisfied")
  }
  n_iter <- posterior_dims[1]
  n_pars <- posterior_dims[2]
  n_sample <- posterior_dims[3]
  
  par_names <- unlist(dimnames(prior_array)[2])
  
  thinner <- seq(from=1, to=n_sample, by=thin)
  ranks <- array(rep(0, n_iter * n_pars), dim=c(n_iter, n_pars))
  dimnames(ranks)[2] <- list(par_names)
  for(i in 1:n_iter){
    for(j in 1:n_pars){
      ranks[i, par_names[j]] <- sum(prior_array[i, par_names[j]] < posterior_array[i, par_names[j], ][thinner])
      print("##############")
      #print(prior_array[i, par_names[j]])
      #print(sum(prior_array[i, par_names[j]] < posterior_array[i, par_names[j], ][thinner]))
    }
  }
  return(ranks)
}

sbc.plot <- function(ranks, par, thin){
  L = dim(ranks)[1] / thin
  CI = qbinom(c(0.005,0.5,0.995), size= dim(ranks)[1], prob = 1/(20))
  print(CI)
  ggplot() + aes(sbc[, par]) + geom_histogram(bins=20) + geom_hline(yintercept = CI, color="black", linetype="dashed")
}