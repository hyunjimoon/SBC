

#' Given prior and posterior samples, generate rank count for each bin.
#'
#' @param prior A named array of dimensions(n_iter, n_pars) where n_iter=number of SBC draws, n_pars the number of parameters. Names should be applied as parameter names.
#' @param posterior A named array of dimensions(n_iter, n_pars, n_samples) where n_iter=number of SBC draws, n_pars the number of parameters, and n_samples the number of samples per SBC draw. Names should be applied as parameter names. Equivalent to posterior.as_draws_array
#' @param thin Integer in which thinning samples will be applied
#' prior array dimensions: (n_iter, n_pars)
#' posterior array dimensions: (n_iter, n_pars, n_sample)
#'
#' @return array of dimensions(n_iter, n_pars)
#' @export
calculate_rank <- function(prior, posterior, thin){

  prior_dims = dim(prior)
  posterior_dims = dim(posterior)

  if (prior_dims[1] != posterior_dims[1] || prior_dims[2] != posterior_dims[2]){
    stop(paste("Dimension mismatch error!",
               "dim(prior)[1] == dim(posterior)[1] && dim(prior)[2] == dim(posterior)[2] must be satisfied",
               paste("prior dimensions:", prior_dims[1], prior_dims[2]),
               paste("posterior dimensions:", posterior_dims[1], posterior_dims[2]),
               "", sep="\n"))
  }
  n_iter <- posterior_dims[1]
  n_pars <- posterior_dims[2]
  n_sample <- posterior_dims[3]

  par_names <- unlist(dimnames(prior)[2])

  thinner <- seq(from=1, to=n_sample, by=thin)
  ranks <- array(rep(0, n_iter * n_pars), dim=c(n_iter, n_pars))
  dimnames(ranks)[2] <- list(par_names)
  for(i in 1:n_iter){
    for(j in 1:n_pars){
      ranks[i, par_names[j]] <- sum(posterior[i, par_names[j], ][thinner] < prior[i, par_names[j]])
    }
  }
  return(ranks)
}
