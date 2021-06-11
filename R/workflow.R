
#' R6 class representing a fluid SBC workflow, boasting a highly customizable pipeline and supports intermediate updates
#' @export
SBCWorkflow <- R6::R6Class("SBCWorkflow",
  public = list(
    #' @field stan_model A CmdStanModel object to use for fitting data
    #' @field sim_function A simulator function, which should return simulated prior and data samples for SBC. Please refer \code{SBCWorkflow$initialize()} for details.
    #' @field calculated_ranks A named list of vectors which contained calculated ranks, or NULL if not calculated.
    #' @field prior_samples A list of named lists which containes simulated prior samples.
    #' @field simulated_y An array with dim(n_iter, n_samples) which contains simulated data.
    #' @field posterior_dimensions A named list containing dimensions of posterior variables. This is used for internal purposes.
    #' @field posterior_samples A list of named lists of posterior samples. names are indexed names, which means multidimensional parameters are decomposed element-wise.
    stan_model = NULL,
    sim_function = NULL,
    calculated_ranks = NULL,
    prior_samples = NULL,  # type is list of named lists. written on simulate
    simulated_y = NULL,  # type is array with dim(n_iter, n_samples). written on simulate
    posterior_dimensions = NULL, # named list containing dimensions of posterior variables. written on fit_model
    posterior_samples = NULL,  # type is list of named lists. written on fit_model

    #' R6 Initializer for SBCWorkflow
    #'
    #' @param stan_model A CmdStanModel object to use for fitting data.
    #' @param sim_function A simulator function. The function must return a
    #'   named list with elements \code{parameters} and \code{generated}.
    #'   Parameters should be a named list containing prior samplers for
    #'   parameters. Generated must be a 1 dimensional vector containing
    #'   simulated y data generated from the \code{parameters} samples.
    #'
    initialize = function(stan_model, sim_function){
      self$stan_model <- stan_model
      self$sim_function <- sim_function
    },

    #' Sample \eqn{\tilde{\theta} \sim P(\theta)} and \eqn{\tilde{y} \sim P(\tilde{y}) | \tilde{\theta} using \code{sim_function}
    #'
    #' @param n Number of prior samples to generate. Equates to number of SBC iterations
    #' @param ... Additional arguments to be passed to \code{sim_function}
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

    #' Sample \eqn{\widehat{\theta} \sim P(\widehat{\theta} | \tilde{y})} using the stan model.
    #'
    #' @param sample_iterations Number of sampling iterations. Equates to number of posterior samples
    #' @param warmup_iterations Number of warmup iteratioins.
    #' @param data List specifying data to be passed to the model. Note that \code{y} will be replaced with simulated \code{y} samples.
    fit_model = function(sample_iterations, warmup_iterations, data=list()){
      if(is.null(self$simulated_y)){
        stop("There are no simulated data available. Please run SBCWorkflow$simulate() first ")
      }
      posterior <- list()  # a list of named lists of posterior draws as r datatypes(vector, array)
      iterations <- dim(self$simulated_y)[1]
      for(n in 1:iterations){
        data[["y"]] <= self$simulated_y[n, ]
        model_fit <- fit_cmdstan_model(self$stan_model, data, sample_iterations, warmup_iterations, 1)
        if(n == 1){
          par_names <- names(self$prior_samples[[1]])
          sample_summary <- as.data.frame(model_fit$summary())
          row.names(sample_summary) <- sample_summary[, "variable"]  # set "variable" column as row names
          par_indexes <- infer_sequential_params(par_names, sample_summary, return_dim_info = TRUE)
          self$posterior_dimensions <- par_indexes[["dims"]]

        }
        par_names <- names(self$prior_samples[[1]])
        drawset <- model_fit$draws(variables=par_names)
        num_indexed_pars <- dim(drawset)[3] # total number of indexed param names
        #print(length(posterior::as_draws_list(drawset)[1]))
        posterior <- append(posterior, posterior::as_draws_list(drawset)[1])
      }
      self$posterior_samples <- posterior
      return(posterior)
    }
  )
)

fit_cmdstan_model <- function(cmdstan_model, data, sampling_iters, warmup_iters, n_chains){
  cmdstan_model$sample(data=data, refresh=0, chains=1, iter_warmup=warmup_iters, iter_sampling=sampling_iters)
}
