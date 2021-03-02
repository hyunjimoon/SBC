infer_sequential_params = function(par_names, summary_data, max_dims=3, return_dim_info=FALSE){
  generate_index_name <- function(name, ...){
    return(paste0(name, "[", paste0(c(...), collapse=","), "]"))  # receive indexes and return stan param name
  }
  return_dims <- list()
  return_names <- list()  # list that holds the final par names, returned

  for(parname in par_names){  # iterate through all received parameter base names
    for(current_dim in max_dims:0){  # start from max dim and descend down to scalar
      if(current_dim == 0){  #  scalar parameter
        return_names <- append(return_names, parname)
        break
      }
      tmp_indexes <- as.integer(rep(1, current_dim))  # initial backtrack attempt index = rep(1, n_dims)

      is_current_dim <- TRUE  # check if value exists for current dimension
      tryCatch(
        {
          mean_ <- summary_data[generate_index_name(parname, tmp_indexes), "mean"]
          stopifnot(!is.na(mean_))  # handle data.frames, which do not raise index errors, but returns NA
        }, error=function(err){is_current_dim <<- FALSE})

      if(isFALSE(is_current_dim)){
        next
      }
      return_names <- append(return_names, generate_index_name(parname, tmp_indexes))  # append first index

      final_indexes <- as.integer(rep(0, current_dim))  # upper bound for each component
      accumulator <- 0  # current component index
      while (sum(tmp_indexes != 0)) {
        skip <- FALSE  # boolean check to handle next requests
        accumulator <- accumulator %% current_dim + 1  # iterate over all components
        if(isTRUE(tmp_indexes[accumulator] == 0)){next}  # skip 0 accumulators
        tryCatch(
          {
            # attempt higher index
            mean_ <- summary_data[generate_index_name(parname, bitwOr(replace(tmp_indexes, accumulator, tmp_indexes[accumulator]+1), final_indexes)), "mean"]
            stopifnot(!is.na(mean_))  # handle data.frames, which do not raise index errors, but returns NA
          },
          error = function(err){  # subscript out of bounds, backtrack and end component
            final_indexes <<- replace(final_indexes, accumulator, tmp_indexes[accumulator])
            tmp_indexes <<- replace(tmp_indexes, accumulator, 0)
            skip <<- TRUE
          }
        )
        if(isTRUE(skip)) {
          next
        }

        tmp_indexes <- replace(tmp_indexes, accumulator, tmp_indexes[accumulator]+1)  # increment backtrack location
        return_names <- append(return_names, generate_index_name(parname, bitwOr(tmp_indexes, final_indexes)))  # tmp_indexes | return_name exists, add to return name
      }
      print(paste(parname, final_indexes))
      return_dims[[parname]] <- final_indexes

      break  # exit current parameter after single iteration
    }
  }
  return(list(names=return_names, dims=return_dims))
}

library(cmdstanr)
library(posterior)

J <- 8
y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)

data = list(J=J, y=y, sigma=sigma)

matrix_test_code <- "
//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[5] eta;
}

transformed parameters {
  matrix[3, 4] mu;
  for(i in 1:3){
    for(j in 1:4){
      mu[i, j] = (i - 1) * 3 + j;
    }
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
}

"

#model <- cmdstan_model(write_stan_file(matrix_test_code))

model = cmdstanr::cmdstan_model("/home/dashadower/git_repos/hyunjimoon/test_SBC/eightschools_ncp.stan")
model_fit <- model$sample(data=data, fixed_param = TRUE)
sample_summary <- as.data.frame(model_fit$summary())  # tibble to data.frame
row.names(sample_summary) <- sample_summary[, "variable"]  # set "variable" column as row names
infer_sequential_params(list("mu", "eta"), sample_summary)



sbc_obj <- SBCModel$new(name="8SNCP", model)

positive_cauchy <- function(...){
  draw <- rcauchy(...)
  return(if(draw >= 0) draw else positive_cauchy(...))
}

hyperpriors <- list(theta_trans=function(){rnorm(J, 0, 1)}, mu=function(){rnorm(1, 0, 5)}, tau=function(){positive_cauchy(1, 0, 5)})

theta_prior <- sbc_obj$sample_theta_tilde(list("theta_trans", "mu", "tau"), 100, hyperpriors)

sampled_y <- sbc_obj$sample_y_tilde(theta_prior, data=data)

sampled_theta <- sbc_obj$sample_theta_bar_y(sampled_y, data=data, pars=list("theta_trans", "tau", "theta"))


res <- model$sample(data=data, iter_warmup=10, iter_sampling=10, chains=1,save_warmup=FALSE, refresh=0, thin=NULL)

#res$draws(variables = list("theta_trans", "tau"))
#extracted <- array(res$draws(variables = list("mu", "tau")), dim = c(10, 2))
