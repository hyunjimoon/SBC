

#' Given prior and posterior samples, generate rank count for each bin.
#'
#' @param prior A named array of dimensions(n_iter, n_pars) where n_iter=number of SBC draws, n_pars the number of parameters. Names should be applied as parameter names.
#' @param posterior A named array of dimensions(n_samples, n_pars, n_iter) where n_iter=number of SBC draws, n_pars the number of parameters, and n_samples the number of samples per SBC draw. Names should be applied as parameter names. Equivalent to posterior.as_draws_array
#' @param thin Integer in which thinning samples will be applied
#' prior array dimensions: (n_iter, n_pars)
#' posterior array dimensions: (n_sample, n_pars, n_iter)
#'
#' @return array of dimensions(n_iter, n_pars)
#' @export
calculate_rank <- function(prior, posterior, thin){

  # Just a hacky conversion from the list representation to make code work right now
  if(is.list(prior)) {
    total_vars <- 0
    for(i in 1:length(prior)) {
      if(length(dim(prior[[i]])) != 2) {
        stop("Multidimensional not supported")
      }
      total_vars <- total_vars + dim(prior[[i]])[2]
    }
    prior_matrix <- matrix(nrow = dim(prior[[1]])[1], ncol = total_vars)
    var_names <- array(NA_character_, total_vars)
    next_index <- 1
    for(i in 1:length(prior)) {
      n_elems <- (dim(prior[[i]])[2])
      if(n_elems == 1) {
        prior_matrix[, next_index] <- prior[[i]][,1]
        var_names[next_index] <- names(prior)[i]
        next_index <- next_index + 1
      } else {
        for(k in 1:n_elems) {
          prior_matrix[, next_index] <- prior[[i]][, k]
          var_names[next_index] <- paste0(names(prior)[i], "[",k,"]")
          next_index <- next_index + 1
        }
      }
    }
    colnames(prior_matrix) <- var_names
    prior <- prior_matrix
  }

  prior_dims = dim(prior)
  posterior_dims = dim(posterior)

  if (prior_dims[1] != posterior_dims[3]){
    stop(paste("Dimension mismatch error!",
               "dim(prior)[1] == dim(posterior)[3] must be satisfied",
               paste("prior dimensions:", prior_dims[1], prior_dims[2]),
               paste("posterior dimensions:", posterior_dims[1], posterior_dims[2], posterior_dims[3]),
               "", sep="\n"))
  }
  n_sample <- posterior_dims[1]
  n_iter <- posterior_dims[3]

  par_names <- intersect(unlist(dimnames(prior)[2]), unlist(dimnames(posterior)[2]))
  n_pars <- length(par_names)
  if(n_pars == 0){
    stop("There isn't a parameter name both present in the prior and posterior column names. SBC cannot continue. length(intersect(prior, posterior)) == 0")
  }

  thinner <- seq(from=1, to=n_sample, by=thin)
  ranks <- array(rep(0, n_iter * n_pars), dim=c(n_iter, n_pars))
  dimnames(ranks)[2] <- list(par_names)
  for(i in 1:n_iter){
    for(j in 1:n_pars){
      ranks[i, par_names[j]] <- sum(posterior[, par_names[j], i][thinner] < prior[i, par_names[j]])
    }
  }
  return(ranks)
}