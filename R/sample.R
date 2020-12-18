
RSTAN_MODEL_CLASS_NAME <- "stanmodel"  # r preprocessor macro equivalent someday?
CMDSTAN_MODEL_CLASS_NAME <- "CmdStanModel"

#' R6 Class Representing a collection of functions to extract data from a single model
#' @export
SBCModel <- R6::R6Class("SBCModel",
  public = list(
    #' @field name Some string to identify your SBCModel
    #' @field stan_model A CmdStanModel or a stanmodel
    #' @field model_type A string containing self$stan_model type. can be either RSTAN_MODEL_CLASS_NAME or CMDSTAN_MODEL_CLASS_NAME
    #' @field hyperpriors A named list containing prior definitions and sampling functions
    name = NULL,
    stan_model = NULL,
    model_type = NULL,
    hyperpriors = NULL,

    #' Initialize a SBCModel with a existing stan Model
    #' @param name Some string to identify your SBCModel
    #' @param stan_model A 'CmdStanModel' or a 'stanmodel' to draw samples
    #' @param hyperpriors A named list containing definitions of prior for parameters. Follows the format:
    #' list(param1=func, param2=func) where func is a function that returns a single value for each parameter
    initialize = function(name, stan_model = NULL, hyperpriors = list()){
      stopifnot(is.character(name))
      stopifnot(CMDSTAN_MODEL_CLASS_NAME %in% class(stan_model) || RSTAN_MODEL_CLASS_NAME %in% class(stan_model))
      stopifnot(length(hyperpriors) > 0)

      if(CMDSTAN_MODEL_CLASS_NAME%in% class(stan_model)){
        self$model_type <- CMDSTAN_MODEL_CLASS_NAME
      }
      else if(RSTAN_MODEL_CLASS_NAME %in% class(stan_model)){
        self$model_type <- RSTAN_MODEL_CLASS_NAME
      }
      else{
        stop(paste("Could not identify stan_model with unknown class type", class(stan_model)))
      }
      self$stan_model <- stan_model
      self$hyperpriors <- hyperpriors
    },

    #' @description
    #' Sample \eqn{\tilde{\theta} ~ P(\theta)} from self$hyperpriors
    #' @param pars List of parameters to draw prior samples.
    #' @param n_iters Integer specifying number of draws.
    #'
    #' @return Array of dimension(n_iters, n_pars) which are sampled prior values.
    sample_theta_tilde = function(pars, n_iters){
      stopifnot(typeof(pars) == "list")

      theta_arr <- array(dim=c(n_iters, length(pars)))
      colnames(theta_arr) <- pars
      for(par_index in 1:length(pars)){
        for(iter_index in 1:n_iters){
          theta_arr[iter_index, pars[[par_index]]] <- self$hyperpriors[[pars[[par_index]]]]()
        }
      }
      return(theta_arr)
    },

    #' @description
    #' Try and use the stan code to sample theta values
    #' We'll set iteration step to 1 and sample from the hyperpriors.
    #'
    #' @param pars_list list of parameter names to draw
    #' @param n_iters integer specifying number of individual draws
    #'
    #' to be implemented(not sure if necessary)
    sample_theta_tilde_stan = function(pars_list, n_iters){
    },

    #' @description
    #' sample \eqn{\tilde{y} ~ p(\tilde{y} | \tilde{\theta})}, or y_tilde ~ P(y_tilde | theta_tilde)
    #' The stan model must specify P(y | theta) within the generated quantities block
    #'
    #' @param theta_arr Array of sampled theta_tilde values (n_iters, n_pars), which is outputted from self$sample_theta_tilde()
    #' @param y_count number of y_tilde samples being generated within generated quantities
    #' @param y_var y_tilde variable name. Will be retrieved from model as so: y_var\[n\] (1<=n<=y_count)
    #' @param data data list to pass to the stan model, in most cases would be dummy data since data will be drawn with fixed params
    #'
    #' @return array of dimension (n_iters, y_count) of sampled y
    sample_y_tilde = function(theta_arr, y_count, y_var="y_rep", data=list()){
      stopifnot(length(dim(theta_arr)) == 2)

      n_iters = dim(theta_arr)[1]
      sample_arr <- array(dim=c(n_iters, y_count))
      for(iter_index in 1:n_iters){
        theta_slice <- as.list(theta_arr[iter_index, ])

        if(self$model_type == CMDSTAN_MODEL_CLASS_NAME){
          # sample for cmdstanr
          model_fit <- self$stan_model$sample(data=data,
                                                 init=list(theta_slice), iter_warmup=1, iter_sampling=1, chains=1,
                                                 parallel_chains=1, save_warmup=FALSE, refresh=0, fixed_param=TRUE)

          sample_summary <- model_fit$summary()
          samples <- lapply(1:y_count, function(x) as.double(sample_summary[sample_summary$variable == paste0(y_var, "[", x, "]"), "mean"]))
          sample_arr[iter_index, ] <- unlist(samples)
        }
        else if(self$model_type == RSTAN_MODEL_CLASS_NAME){
          # sample for rstan
          model_fit <- rstan::sampling(self$stan_model, data=data,
                                       pars=as.vector(c(dimnames(theta_arr)[[2]], y_var)), chains=1, algorithm="Fixed_param",
                                       init=list(theta_slice), cores=1, show_messages=FALSE, refresh=0)

          sample_summary <- rstan::summary(model_fit)$summary
          samples <- lapply(1:y_count, function(x) sample_summary[paste0(y_var, "[", x, "]"), "mean"])
          sample_arr[iter_index, ] <- unlist(samples)
        }
      }
      rm(model_fit, sample_summary, samples)  # cleanup
      return(sample_arr)
    },

    #' @description
    #' sample \eqn{\theta ~ p(\theta | \tilde{y})}
    #' The stan model must have \eqn{\theta} defined as parameters
    #'
    #' @param y_sample_arr array of dimension (n_iters, y_count) of sampled y. Same as output of self$sample_y_tilde()
    #' @param data list of additional data to pass to the stan model. Note that "y" will be overwritten with draws from y_sample_arr.
    #' @param pars list of parameters of interest.
    #' @param fit_iter number of model iterations.
    #'
    #' @return: array of dimension (n_iters, n_pars, fit_iter) of posterior parameter draws
    sample_theta_bar_y = function(y_sample_arr, data=list(), pars=list(), fit_iter=200){
      stopifnot(length(pars) > 0)

      n_iters = dim(y_sample_arr)[1]
      draw_arr <- array(dim=c(n_iters, length(pars), fit_iter))
      if(typeof(pars) == "list"){
        pars <- unlist(pars)
      }

      for(iter_index in 1:n_iters){
        data[["y"]] <- y_sample_arr[iter_index, ]  # insert "y" with sample slice vector

        if(self$model_type == CMDSTAN_MODEL_CLASS_NAME){
          model_fit <- self$stan_model$sample(data=data, iter_warmup=fit_iter, iter_sampling=fit_iter, chains=1,
                                              save_warmup=FALSE, refresh=0, thin=NULL)

          draw_arr[iter_index, , ] <- aperm(model_fit$draws(variables=pars)[, 1, ])  # arrays are filled row first, so we permutate once
          # https://github.com/stan-dev/cmdstanr/issues/346
          # currently throws an error if you try to retrieve only a single numeric variable, just add lp until fix is released
        }

        else if(self$model_type == RSTAN_MODEL_CLASS_NAME){
          model_fit <- rstan::sampling(self$stan_model, data=data, pars=pars, chains=1, iter=fit_iter*2, warmup=fit_iter,
                                       refresh=0)

          draw_arr[iter_index, , ] <- aperm(posterior::as_draws_array(rstan::extract(model_fit, pars=pars, permuted=FALSE, inc_warmup=FALSE)))
        }
      }
      if(dim(draw_arr)[2] != length(pars)){
        warning("Couldn't rename array dimnames with corresponding parameters.
                This probably means a sequential parameter is included, and SBC can't verify the identities of parameter samples.
                If you would like to use SBC plotting features, you need to manually define dimnames for each parameter in dimension 2")
      }
      else{
        dimnames(draw_arr)[[2]] <- pars  # TODO: handle sequential parameters
      }
      rm(model_fit)
      return(draw_arr)
    },

    #' @description
    #' sample \eqn{\tilde{y}^{*} ~ Bootstrap(\tilde{y})}
    #' given a single vector of data y_tilde, bootstrap sample multiple values of y_tilde_star\[i\] from bootstrap(y_tilde)
    #'
    #' @param y_sample_vector vector of y samples.
    #' @param n_iters number of bootstrap samples to draw
    #'
    #' @return array of dimension (n_iters, y_count) of bootstrap sampled y
    sample_bootstrap_y_tilde = function(y_sample_vector, n_iters=2000){
      sample_cnt = length(y_sample_vector)
      y_star_arr <- array(dim=c(n_iters, sample_cnt))
      for(iter_index in 1:n_iters){
        y_star_arr[iter_index, ] <- sample(y_sample_vector, replace=TRUE)  # sample from y_tilde
      }
      return(y_star_arr)
    }
  )
)
