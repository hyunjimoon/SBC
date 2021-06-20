
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
    #' @field thin_factor Non-negative integer representing posterior thinning interval
    stan_model = NULL,
    sim_function = NULL,
    calculated_ranks = NULL,
    prior_samples = NULL,  # type is list of named lists. written on simulate
    simulated_y = NULL,  # type is array with dim(n_iter, n_samples). written on simulate
    posterior_dimensions = NULL, # named list containing dimensions of posterior variables. written on fit_model
    posterior_samples = NULL,  # type is list of named lists. written on fit_model
    thin_factor = NULL, # type is integer. written on initialize and fit_model

    #' R6 Initializer for SBCWorkflow
    #'
    #' @param stan_model A CmdStanModel object to use for fitting data.
    #' @param sim_function A simulator function. The function must return a
    #'   named list with elements \code{parameters} and \code{generated}.
    #'   Parameters should be a named list containing prior samplers for
    #'   parameters. Generated must be a 1 dimensional vector containing
    #'   simulated y data generated from the \code{parameters} samples.
    #'
    #' @return None
    initialize = function(stan_model, sim_function, thin_factor=3){
      self$stan_model <- stan_model
      self$sim_function <- sim_function
      self$thin_factor <- thin_factor
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
    #' @param thin_factor Non-negative integer in which thinning should be applied. 1 equates to no thinning being done.
    #'
    #' @return A list of named lists containing posterior samples for each parameter.
    fit_model = function(sample_iterations, warmup_iterations, data=list(), thin_factor=3){
      self$thin_factor <- thin_factor
      if(is.null(self$simulated_y)){
        stop("There are no simulated data available. Please run SBCWorkflow$simulate() first ")
      }
      posterior <- NULL  # a list of named lists of posterior draws as r datatypes(vector, array)
      iterations <- dim(self$simulated_y)[1]
      for(n in 1:iterations){
        data[["y"]] <= self$simulated_y[n, ]
        model_fit <- fit_cmdstan_model(self$stan_model, data, sample_iterations * thin_factor, warmup_iterations , 1, thin=thin_factor)
        par_names <- names(self$prior_samples[[1]])
        drawset <- model_fit$draws(variables=par_names)
        if(n == 1){
          sample_summary <- as.data.frame(model_fit$summary())
          row.names(sample_summary) <- sample_summary[, "variable"]  # set "variable" column as row names
          par_indexes <- infer_sequential_params(par_names, sample_summary, return_dim_info = TRUE)
          self$posterior_dimensions <- par_indexes[["dims"]]
          posterior <- list(posterior::as_draws_list(drawset)[[1]])  # re-initialize array
        }
        num_indexed_pars <- dim(drawset)[3] # total number of indexed param names
        posterior <- append(posterior, list(posterior::as_draws_list(drawset)[[1]]))
      }
      self$posterior_samples <- posterior
      return(posterior)
    },

    #' Calculate rank statistics for a given parameter name
    #'
    #' @param param list of parameter names to calculate. If not given, calculate for all available parameters.
    #'
    calculate_rank = function(param=NULL){
      if(is.null(self$posterior_samples)){
        stop("There are no posterior samples available. Please run SBCWorkflow$fit_model() first.")
      }
      if(is.null(param)){
        param <- names(self$posterior_samples[[1]])
      }

      n_iters <- length(self$prior_samples)
      n_params <- length(param)
      decomposed_prior <- decompose_structure_to_par_list(self$prior_samples[[1]])
      base_param_names <- unlist(lapply(param, function(n){if(grepl("\\[", n)) sub("\\[.*", "", n) else n}))  # extract base parameter names(A[1, 2] -> A)
      if(length(decomposed_prior) != n_params){
        stop("Number of prior and posterior parameters are mismatched. Please ensure each parameters' dimensions are corretly defined.")
      }

      ranks <- array(rep(0, n_iters * n_params), dim=c(n_iters, n_params))
      dimnames(ranks)[2] <- list(names(decomposed_prior))

      for(i in 1:n_iters){
        prior <- self$prior_samples[[i]]
        decomposed_prior <- decompose_structure_to_par_list(prior)

        for(j in 1:n_params){
          regex_res <- gregexpr("(?<=\\[)(.+)(?=\\])", param[[j]], perl = TRUE)[[1]]
          capture_start <- attr(regex_res, "capture.start")
          capture_length <- attr(regex_res, "capture.length")
          if(capture_start != -1){
            index <- as.integer(unlist(substr(param[[j]], capture_start, capture_start + capture_length - 1)))
            ranks[i, param[[j]]] <- sum(self$posterior_samples[[i]][[param[[j]]]] < access_element_by_index(self$prior_samples[[i]][[base_param_names[[j]]]], index))
          }
          else{
            t1 <- self$posterior_samples[[i]][[param[[j]]]]
            t2 <- self$prior_samples[[i]][[param[[j]]]]
            #ranks[i, param[[j]]] <- sum(self$posterior_samples[[i]][[param[[j]]]] < self$prior_samples[[i]][[param[[j]]]])
            ranks[i, param[[j]]] <- sum(t1 < t2)
          }
        }
      }
      return(ranks)
    }
  )
)

fit_cmdstan_model <- function(cmdstan_model, data, sampling_iters, warmup_iters, n_chains, ...){
  return(do.call(cmdstan_model$sample, list(data=data, iter_sampling=sampling_iters, iter_warmup=warmup_iters, chains=n_chains, ...)))
  #cmdstan_model$sample(data=data, refresh=0, chains=1, iter_warmup=warmup_iters, iter_sampling=sampling_iters)
}
