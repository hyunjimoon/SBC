library(R6)

RSTAN_MODEL_CLASS_NAME <- "stanmodel"  # sucks that R doesn't have 
CMDSTAN_MODEL_CLASS_NAME <- "CmdStanModel"

SBCModel <- R6Class("SBCModel", list(
  # @field name: Some string to identify your SBCModel
  # @field stan_model: A CmdStanModel or a stanmodel
  # @field model_type: A string containing self$stan_model type
  # @field hyperpriors: A named list containing prior definitions and sampling functions
  name = NULL,
  stan_model = NULL,
  model_type = NULL,
  hyperpriors = NULL,
  initialize = function(name, stan_model = NULL, hyperpriors = list()){
    # @param name: Some string to identify your SBCModel
    # @param stan_model: A 'CmdStanModel' or a 'stanmodel' to draw samples
    # @param hyperpriors: A named list containing definitions of prior for parameters. Follows the format:
    # list(param1=func, param2=func) where func is a function that returns a single value for each parameter
    stopifnot(is.character(name))
    stopifnot(CMDSTAN_MODEL_CLASS_NAME %in% class(stan_model) || RSTAN_MODEL_CLASS_NAME %in% class(stan_model))
    stopifnot(length(hyperpriors) > 0)

    if(CMDSTAN_MODEL_CLASS_NAME%in% class(stan_model)){
      self$model_type <- CMDSTAN_MODEL_CLASS_NAME
    }
    else if(RSTAN_MODEL_CLASS_NAME %in% class(stan_model)){
      self$model_type <- RSTAN_MODEL_CLASS_NAME
    }
    self$stan_model <- stan_model
    self$hyperpriors <- hyperpriors
  },
  sample_theta_tilde = function(pars, iters){
    # @param pars: List of parameters to draw prior samples.
    # @param iters: Integer specifying number of draws.
    #
    # returns: Array of dimension(n_iters, n_pars) which are sampled prior values.
    stopifnot(typeof(pars) == "list")

    theta_arr <- array(dim=c(iters, length(pars)))
    colnames(theta_arr) <- pars
    for(par_index in 1:length(pars)){
      for(iter_index in 1:iters){
        theta_arr[iter_index, pars[[par_index]]] <- self$hyperpriors[[pars[[par_index]]]]()
      }
    }
    return(theta_arr)
  },
  sample_y_tilde = function(theta_arr, y_count, y_var="y_rep", data=list()){
    # sample $\tilde{y} ~ p(\tilde{y} | \tilde{\theta})$, or y_tilde ~ P(y_tilde | theta_tilde)
    # The stan model must specify P(y | theta) within the generated quantities block
    #
    # @param theta_arr: Array of sampled theta_tilde values, which is outputted from self$sample_theta_tilde()
    # @param y_count: number of y_tilde samples being generated within generated quantities
    # @param y_var: y_tilde variable name. Will be retrieved from model as so: y_var[n] (1<=n<=y_count)
    # @param data: data list to pass to the stan model, in most cases is irrelevant since sampling will be done with fixed params
    #
    # returns: array of length (n_iters, y_count) of sampled y
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
  approx_theta_bar_y = function(y_sample_arr){
  }
))