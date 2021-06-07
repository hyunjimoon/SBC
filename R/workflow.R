
#' R6 class representing a fluid SBC workflow, boasting a highly customizable pipeline and supports intermediate updates
#' @export
SBCWorkflow <- R6::R6Class("SBCWorkflow",
  public = list(
    #' @field name Some string that identifies this workflow (for your convenience)
    name = NULL,
    stan_model = NULL,
    sim_function = NULL,
    calculated_ranks = NULL,
    prior_samples = NULL,  # type is list of named lists
    simulated_y = NULL,  # type is array

    initialize = function(stan_model, sim_function, ...){
      self$stan_model <- stan_model
      self$sim_function <- sim_function
    },

    simulate = function(n, ...){  # number of simulation draws should be at least 1000
      prior_list <- list()
      sim_y <- NA
      for(i in 1:n){
        sim_list <- do.call(self$sim_function, list(...))  # returned type is a list
        priors <- sim_list[["parameters"]]
        prior_list <- append(prior_list, list(priors))
        sample_y <- sim_list[["generated"]]
        if(i == 1){
          sim_y <- array(rep(NA, n * length(sample_y)), c(n, length(sample_y)))
        }
        sim_y[i, ] <- unlist(sample_y)
      }
      self$prior_samples <- prior_list
      self$simulated_y <- sim_y

    },

    fit_model = function(prior_samples, sample_iterations, warmup, data=list()){
      fit_cmdstan_model(self$stan_model, data, sample_iterations, warmup, 1)
    },

    sample_data = function(priors, ...){
      sim_y <- array()
      for(i in 1:length(priors)){
        sim_row <- do.call(self$data_sim_function, list(priors[[i]], list(...)))
        if(i == 1){
          sim_y <<- array(rep(NA, length(priors) * length(sim_row)), c(length(priors), length(sim_row)))
        }
        sim_y[i, ] <<- unlist(sim_row)
      }
      return(sim_y)
    }

  ))

fit_cmdstan_model <- function(cmdstan_model, data, sampling_iters, warmup_iters, n_chains){
  cmdstan_model$sample(data=data, refresh=0, chains=1, iter_warmup=warmup_iters, iter_sampling=sampling_iters)
}


convert_cmdstanfit <- function(cmdstan_fit){
  sample_summary <- as.data.frame(cmdstan_fit$summary())
}

multidim_index_to_contiguous <- function(index, dims){
  arrayInd(index, dim)
}
