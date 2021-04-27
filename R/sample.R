
RSTAN_MODEL_CLASS_NAME <- "stanmodel"  # r preprocessor macro equivalent someday?
CMDSTAN_MODEL_CLASS_NAME <- "CmdStanModel"

#' R6 Class Representing a collection of functions to extract data from a single model
#' @export
SBCModel <- R6::R6Class("SBCModel",
  public = list(
    #' @field name Some string to identify your SBCModel
    #' @field stan_model A CmdStanModel or a stanmodel
    #' @field model_type A string containing self$stan_model type. can be either RSTAN_MODEL_CLASS_NAME or CMDSTAN_MODEL_CLASS_NAME
    #' @field parameter_dims A list containing dimension info for parameters. Parameter names are original, unsuffixed names.
    #' @field sum_suffix Characters representing the suffix added to simulated parameters
    name = NULL,
    stan_model = NULL,
    model_type = NULL,
    parameter_dims = NULL,
    sim_suffix = NULL,


    #' Initialize a SBCModel with a existing stan Model
    #' @param name Some string to identify your SBCModel
    #' @param stan_model A 'CmdStanModel' or a 'stanmodel' to draw samples
    initialize = function(name, stan_model = NULL){
      stopifnot(is.character(name))
      stopifnot(CMDSTAN_MODEL_CLASS_NAME %in% class(stan_model) || RSTAN_MODEL_CLASS_NAME %in% class(stan_model))

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
    },

    #' @description
    #' Sample \eqn{\tilde{\theta} ~ P(\theta)} (Prior predictive sampling)
    #' @param pars List of parameters to draw prior samples.
    #' @param n_iters Integer specifying number of draws.
    #' @param hyperpriors A named list containing prior definitions and sampling functions. Follows the format:
    #' list(param1=func, param2=func) where func is a function that returns a single value for each parameter
    #'
    #' @return List of length(pars) where each element is an array of dimension(n_iters, length(par\[i\])) containing sampled prior values.
    sample_theta_tilde = function(pars, n_iters, hyperpriors=list()){
      stopifnot(typeof(pars) == "list")
      stopifnot(length(hyperpriors) > 0)

      theta_list <- list()
      # First iteration: determine the dimension of each parameter and initialize arrays
      for(par_index in 1:length(pars)){
        sampled_par <- hyperpriors[[pars[[par_index]]]]()
        par_dims <- length(sampled_par)
        theta_list[[pars[[par_index]]]] <- array(dim=c(n_iters, par_dims))
        theta_list[[pars[[par_index]]]][1, ] <- sampled_par
      }

      # Second iteration and onwards: fill up the array with samples
      if (n_iters > 1){
        for(par_index in 1:length(pars)){
          for(iter_index in 2:n_iters){
            theta_list[[pars[[par_index]]]][iter_index, ] <- hyperpriors[[pars[[par_index]]]]()
          }
        }
      }
      return(theta_list)
    },

    #' @description
    #' Try and use the stan code to sample theta values(prior predictive sampling).
    #' This is useful if you have prior simulation defined in generated quantities.
    #' Note that this does *not* use the model, since Stan doesn't directly sample from the prior.
    #'
    #' @param pars_list list of parameter names to draw
    #' @param suffix additional suffix added to simulated parameter names. default is "_". For example,
    #' if the parameter is called "mu", "mu_" will be sampled as the prior predictive distribution.
    #' @param n_iters integer specifying number of individual draws
    #' @param data additional data for the model, if necessary
    #'
    #' @return List of length(pars) where each element is an array of dimension(n_iters, length(par\[i\])) containing sampled prior values.
    sample_theta_tilde_stan = function(pars_list, n_iters, data=list(), suffix="_"){
      ## TODO: Rewrite to match return type of sample_theta_tilde (list instead of array)
      stopifnot(length(pars_list) > 0)
      self$sim_suffix <- suffix
      suffixed_pars_list = lapply(pars_list, function(x){paste0(x, suffix)})  # actual parameter names
      return_arr <- NULL
      for(iter_count in 1:n_iters){
        if(self$model_type == CMDSTAN_MODEL_CLASS_NAME){
          # explore parameters for cmdstan
          model_fit <- self$stan_model$sample(iter_warmup=1, iter_sampling=1, chains=1, fixed_param=TRUE,
                                              data=data, parallel_chains=1)

          sample_summary <- as.data.frame(model_fit$summary())  # tibble to data.frame
          row.names(sample_summary) <- sample_summary[, "variable"]  # set "variable" column as row names

        }
        else if(self$model_Type == RSTAN_MODEL_CLASS_NAME){
          # explore parameters for rstan
          model_fit <- rstan::sampling(self$stan_model, data=data,
                                       chains=1, algorithm="Fixed_param",
                                       iter=2, warmup=1, cores=1, show_messages=FALSE, refresh=0)

          sample_summary <- rstan::summary(model_fit)$summary
        }

        suffix_indexed_param_names <- self$infer_sequential_params(suffixed_pars_list, sample_summary, return_dim_info = TRUE)

        self$parameter_dims <- suffix_indexed_param_names[["dims"]]
        names(self$parameter_dims) <- lapply(names(self$parameter_dims), function(n){stringi::stri_replace_last_fixed(n, suffix, "")})

        suffix_indexed_param_names <- suffix_indexed_param_names[["names"]]
        indexed_param_names <- lapply(suffix_indexed_param_names, function(n){stringi::stri_replace_last_fixed(n, suffix, "")})
        if(is.null(return_arr)){  # initialize array
          return_arr <- array(dim=c(n_iters, length(suffix_indexed_param_names)))
          colnames(return_arr) <- indexed_param_names
          lapply(suffix_indexed_param_names, function(col_name){  # need to run assignment
            return_arr[iter_count, stringi::stri_replace_last_fixed(col_name, suffix, "")] <<- sample_summary[col_name, "mean"]  # retrieve mean
          })
        }
        else{  # fill values within column
          lapply(suffix_indexed_param_names, function(col_name){
            return_arr[iter_count, stringi::stri_replace_last_fixed(col_name, suffix, "")] <<- sample_summary[col_name, "mean"]  # retrieve mean
          })
        }
      }
      # any(unlist(lapply(theta_list, function(n){any(is.na(n))})))
      if(!all(!is.na(return_arr))){
        warning("Some samples contain NA as values. This means something went terribly wrong during sampling, or SBC is
                failing to identify one of your prior parameters' names. Please check you have parameters with suffixes defined within
                your model.")
      }
      return(return_arr)
    },

    #' @description
    #' sample \eqn{\tilde{y} ~ p(\tilde{y} | \tilde{\theta})}, or y_tilde ~ P(y_tilde | theta_tilde),
    #' which is the Posterior Predictive Distribution.
    #' The stan model must specify P(y | theta) within the generated quantities block
    #'
    #' @param theta_arr Array of sampled prior theta values, which is outputted from self$sample_theta_tilde_stan()
    #' @param y_var y_tilde variable name. Will be retrieved from model as so: y_\[n\] (1<=n<=y_count)
    #' @param data data list to pass to the stan model, in most cases would be dummy data since data will be drawn with fixed params
    #'
    #' @return array of dimension (n_iters, y_dim, ...) of sampled y. Each row is a sample vector for a single parameter vector
    sample_y_tilde = function(theta_arr, y_var="y_", data=list()){

      #n_iters = dim_t(theta_list[[names(theta_list)[1]]])[1]
      n_iters = dim(theta_arr)[[1]]
      for(iter_index in 1:n_iters){
        theta_slice <- as.list(theta_arr[iter_index, ])
        #theta_slice <- lapply(theta_list, function(x){x[iter_index, ]})
        if(self$model_type == CMDSTAN_MODEL_CLASS_NAME){
          # sample for cmdstanr
          model_fit <- self$stan_model$sample(data=data,
                                                 init=list(theta_slice), iter_warmup=1, iter_sampling=1, chains=1,
                                                 parallel_chains=1, save_warmup=FALSE, refresh=0, fixed_param=TRUE)

          sample_summary <- as.data.frame(model_fit$summary())  # tibble to data.frame
          row.names(sample_summary) <- sample_summary[, "variable"]  # set "variable" column as row names
          y_indexes = self$infer_sequential_params(list(y_var), sample_summary, return_dim_info = TRUE)
          if(iter_index == 1){
            sample_arr <<- array(dim=c(n_iters, y_indexes[["dims"]][[y_var]]))
          }
          samples <- unlist(lapply(y_indexes[["names"]], function(x) {sample_summary[x, "mean"]}))
          sample_arr[iter_index, ] <- samples
        }
        else if(self$model_type == RSTAN_MODEL_CLASS_NAME){
          # sample for rstan
          model_fit <- rstan::sampling(self$stan_model, data=data,
                                       pars=as.vector(c(names(theta_list), y_var)), chains=1, algorithm="Fixed_param",
                                       init=list(theta_slice), iter=2, warmup=1, cores=1, show_messages=FALSE, refresh=0)

          sample_summary <- rstan::summary(model_fit)$summary
          y_indexes = self$infer_sequential_params(list(y_var), sample_summary, return_dim_info = TRUE)
          if(iter_index == 1){
            sample_arr <<- array(dim=c(n_iters, y_indexes[["dims"]][[y_var]]))
          }
          samples <- unlist(lapply(y_indexes["names"], function(x) {as.double(sample_summary[x, "mean"])}))
          sample_arr[iter_index, ] <- samples
        }
      }
      rm(model_fit, sample_summary, samples)  # cleanup
      if(!all(!is.na(sample_arr))){
        warning(paste0("NA values are present within extracted samples. This means something went terribly wrong during sampling, or SBC is
                failing to find the simulated data variable. Please check that you have the generated data variable, which is set to be defined as (", y_var, ")."))
      }
      self$par_list_to_structure(theta_slice)
      return(sample_arr)
    },

    #' @description
    #' Sample \eqn{\theta ~ p(\theta | \tilde{y})}, to retrieve parameters from samples.
    #' The stan model must have \eqn{\theta} defined as parameters
    #'
    #' @param y_sample_arr array of dimension (n_iters, y_dim, ...) of sampled y. Same as output of self$sample_y_tilde()
    #' @param data list of additional data to pass to the stan model. Note that "y" will be overwritten with draws from y_sample_arr.
    #' @param pars list of parameters of interest.
    #' @param fit_iter number of model iterations, which equates to the number of posterior draws per sample set.
    #'
    #' @return: array of dimension (fit_iter, n_pars, n_iters) of posterior parameter draws
    sample_theta_bar_y = function(y_sample_arr, data=list(), pars=list(), fit_iter=200){
      stopifnot(length(pars) > 0)

      n_iters = dim(y_sample_arr)[1]
      if(typeof(pars) == "list"){
        pars <- unlist(pars)
      }

      draw_arr <- NULL
      num_indexed_pars <- NULL

      # first iteration - determine length of parameter vector and dimnames
      data[["y"]] <- y_sample_arr[1, ]  # insert "y" with sample slice vector

      if(self$model_type == CMDSTAN_MODEL_CLASS_NAME){
        model_fit <- self$stan_model$sample(data=data, iter_warmup=fit_iter, iter_sampling=fit_iter, chains=1,
                                            save_warmup=FALSE, refresh=0, thin=NULL)

        par_arr <- model_fit$draws(variables=pars)
        num_indexed_pars <- dim(par_arr)[3]
        draw_arr <- array(dim=c(fit_iter, num_indexed_pars, n_iters), dimnames = list(c(1:fit_iter), dimnames(par_arr)[["variable"]], c(1:n_iters)))
        draw_arr[, , 1] <- array(par_arr, dim=c(fit_iter, num_indexed_pars))
      }

      else{
        # rstan
        # TODO: UPDATE AND TEST FOR RSTAN
      }

      # onwards: fill up the array
      for(iter_index in 2:n_iters){
        data[["y"]] <- y_sample_arr[iter_index, ]  # insert "y" with sample slice vector

        if(self$model_type == CMDSTAN_MODEL_CLASS_NAME){
          model_fit <- self$stan_model$sample(data=data, iter_warmup=fit_iter, iter_sampling=fit_iter, chains=1,
                                              save_warmup=FALSE, refresh=0, thin=NULL)

          draw_arr[, , iter_index] <- array(model_fit$draws(variables=pars), dim=c(fit_iter, num_indexed_pars))
        }

        else if(self$model_type == RSTAN_MODEL_CLASS_NAME){
          # TODO: UPDATE AND TEST FOR RSTAN
          model_fit <- rstan::sampling(self$stan_model, data=data, pars=pars, chains=1, iter=fit_iter*2, warmup=fit_iter,
                                       refresh=0)

          draw_arr[, , iter_index] <- posterior::as_draws_array(rstan::extract(model_fit, pars=pars, permuted=FALSE, inc_warmup=FALSE))
        }
      }
      if(dim(draw_arr)[2] != num_indexed_pars){
        warning("Couldn't rename array dimnames with corresponding parameters.
                This probably means a sequential parameter is included, and SBC can't verify the identities of parameter samples.
                If you would like to use SBC plotting features, you need to manually define dimnames for each parameter in dimension 2")
      }
      #else{
        #dimnames(draw_arr)[[2]] <- pars  #TODO: handle sequential parameters
      #}
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
    },

    #' @description
    #' High level function that automatically runs all low level sampling functions, and returns just posterior rank statistics.
    #' Output is same as SBC::calculte_rank()
    #'
    #' @param priors named list of functions which define how each parameter's prior is sampled. If empty, will be
    #' sampled from stan with the assumption that it's defined within generated quantities.
    #' @param pars list of base parameter names of interest(no indexes)
    #' @param prior_suffix additional suffix added to simulated parameter names. default is "_"
    #' @param n_iters number of SBC iterations to run. default is 20
    #' @param n_fits number of iterations to fit theta. default is 200
    #' @param data list of additional data necessary for the model to run
    #' @param y_gq_name name of posterior predictive variable for y. default is "y_"
    #' @param thin sample thinning interval. default is 3
    #'
    #' @return array of dimension (n_iters, n_pars) of posterior rank statistics
    sample_all = function(priors=list(), pars=list(), prior_suffix="_", n_iters=20, n_fits=200, data=list(), y_gq_name="y_", thin=3){
      stopifnot(length(pars) > 0)

      if(length(priors) == 0){
        theta_sample <- self$sample_theta_tilde_stan(pars, n_iters, data, prior_suffix)
      }
      else{
        theta_sample <- self$sample_theta_tilde(pars, n_iters)
      }

      sampled_y <- self$sample_y_tilde(theta_sample, data=data, y_var=y_gq_name)
      theta_post <- self$sample_theta_bar_y(sampled_y, data=data, pars=pars, fit_iter=n_fits * thin)
      return(SBC::calculate_rank(theta_sample, theta_post, thin))
    },

    #' @description
    #' Given list of base parameter names and stan summary matrix, return list of all individual index parameter names
    #' For example, if the model contained a vector\[5\] parname as a parameter, you can run infer_sequential_param("parname", ...)
    #' and it will return a list of strings "parname\[1\]", "parname\[2\]", ...
    #' This function will also search for multidimensional array types, up to dimension (max_dims).
    #'
    #' @param par_names list of base parameter names to find indexes
    #' @param summary_data stansummary matrix that has row names as indexed parameter names, and at least "mean" in its column. It
    #' must raise a subscript error if a row that doesn't exist is indexed. Same as output of rstan::summary
    #' @param max_dims Maximum dimension to search
    #' @param return_dim_info Boolean indicating whether to additionally return dimension info. Default is FALSE
    #'
    #' @return list of strings, that represent indexed parameter names, or a list containing string and dimensions for each parameter if return_dim_info is TRUE
    infer_sequential_params = function(par_names, summary_data, max_dims=3, return_dim_info=FALSE){
      generate_index_name <- function(name, ...){
        return(paste0(name, "[", paste0(c(...), collapse=","), "]"))  # receive indexes and return stan param name
      }
      return_dims <- list()  # list that holds the dimensions for each parameter
      return_names <- list()  # list that holds the final par names, returned

      for(parname in par_names){  # iterate through all received parameter base names
        for(current_dim in max_dims:0){  # start from max_dims and descend down to scalar
          if(current_dim == 0){  #  scalar parameter
            return_names <- append(return_names, parname)
            return_dims[[parname]] <- c(1)
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

          while(sum(tmp_indexes) != 0) {
            skip <- FALSE  # boolean check to handle next requests
            accumulator <- accumulator %% current_dim + 1  # iterate over all components
            if(isTRUE(tmp_indexes[accumulator] == 0)){next}  # skip 0 accumulators
            tryCatch(
              {
                # attempt higher index
                print(paste("Attempting", generate_index_name(parname, bitwOr(replace(tmp_indexes, accumulator, tmp_indexes[accumulator]), final_indexes))))
                #mean_ <- summary_data[generate_index_name(parname, bitwOr(replace(tmp_indexes, accumulator, tmp_indexes[accumulator]+1), final_indexes)), "mean"]
                mean_ <- summary_data[generate_index_name(parname, bitwOr(replace(tmp_indexes, accumulator, tmp_indexes[accumulator]), final_indexes)), "mean"]
                stopifnot(!is.na(mean_))  # handle data.frames, which do not raise index errors, but returns NA
              },
              error = function(err){  # subscript out of bounds, backtrack and end component
                print(paste("Failed", generate_index_name(parname, bitwOr(replace(tmp_indexes, accumulator, tmp_indexes[accumulator]), final_indexes))))
                if(current_dim == 1){
                  final_indexes <<- replace(final_indexes, accumulator, tmp_indexes[accumulator] - 1)
                }
                else{
                  final_indexes <<- replace(final_indexes, accumulator, tmp_indexes[accumulator])
                }
                tmp_indexes <<- replace(tmp_indexes, accumulator, 0)
                skip <<- TRUE
              }
            )
            if(isTRUE(skip)) {
              next
            }

            tmp_indexes <- replace(tmp_indexes, accumulator, tmp_indexes[accumulator]+1)  # increment backtrack location
            #return_names <- append(return_names, generate_index_name(parname, bitwOr(tmp_indexes, final_indexes)))  # tmp_indexes | return_name exists, add to return name
          }

          final_indexes <- rev(final_indexes)
          return_dims[[parname]] <- final_indexes
          for(i in 2:prod(final_indexes)){
            return_names <- append(return_names, generate_index_name(parname, arrayInd(i, final_indexes)))
          }
          break  # exit current parameter after single iteration
        }
      }
      if(return_dim_info){
        return(list(names=return_names, dims=return_dims))
      }
      else{
        return(return_names)
      }
    },

    #' @description
    #' Given a named list, where indexes are parameter names with bracketed indexes(in the form of mu\[1\], mu\[2\]),
    #' convert indexed parameters into actual R vectors and matrices.
    #'
    #' @param named_list A named list containing bracketed indexed parameter names as indexes.
    #'
    #' @return A list with bracketed names removed, and with adequate vector/matrix structures.
    par_list_to_structure = function(named_list){
      unique_par_names <- unlist(unique(lapply(names(named_list), function(n){if(grepl("\\[", n)) sub("\\[.*", "", n) else n})))  # extract base parameter names(A[1, 2] -> A)
      generate_index_name <- function(name, ...){
        return(paste0(name, "[", paste0(c(...), collapse=","), "]"))  # receive indexes and return stan param name
      }
      return_list = list()
      for(i in 1:length(unique_par_names)){
        suffixed <- FALSE
        print("-------------------------")
        suffixed_parname <- unique_par_names[[i]]
        base_parname <- suffixed_parname
        if(stringi::stri_sub(base_parname, -length(self$sim_suffix), -1) == self$sim_suffix){
          # If suffixed, remove suffix from base_parname
          suffixed <- TRUE
          base_parname <- stringi::stri_replace_last_fixed(base_parname, self$sim_suffix, "")
        }

        dim_length <- length(self$parameter_dims[[base_parname]])
        par_dims <- self$parameter_dims[[base_parname]]
        print(paste("base_parname", base_parname))
        print(paste("dim_length", dim_length))
        print(paste("par_dims", par_dims))


        if(dim_length == 1){
          if(par_dims == 1){  # scalar parameter
            return_list[[if(suffixed) suffixed_parname else base_parname]] <-named_list[[base_parname]]
          }
          else{  # vector or 1d array
            return_list[[if(suffixed) suffixed_parname else base_parname]] <- rep(0, self$parameter_dims[[i]])
            for(comp_index in 1:(self$parameter_dims[[i]])){
              return_list[[if(suffixed) suffixed_parname else base_parname]][[comp_index]] <- named_list[[generate_index_name(base_parname, c(comp_index))]]
            }
          }
        }
        else if(dim_length > 1){
          print("entry")
          n_components <- prod(par_dims)
          print(paste("n_components", n_components))
          linear_array <- rep(NA, n_components)
          for(i in 1:n_components){
            arr_ind <- arrayInd(i, par_dims)[1, ]
            print(paste("arr_ind", arr_ind))
            print(paste("index:", generate_index_name(base_parname, arr_ind), "actual_value:", as.numeric(named_list[[generate_index_name(base_parname, arr_ind)]])))
            linear_array[[i]] <- as.numeric(named_list[[generate_index_name(base_parname, arr_ind)]])
          }
          return_list[[if(suffixed) suffixed_parname else base_parname]] <- array(as.numeric(unlist(linear_array)), dim=par_dims)
        }
        print(return_list[[if(suffixed) suffixed_parname else base_parname]])
      }
      print("all done!")
      return(return_list)
    }
  )
)


##########################
# some helper functions

dim_t <- function(data_t){
  # "templated" dim() function
  if(is.array(data_t) | is.matrix(data_t) | is.data.frame(data_t)){
    return(dim(data_t))
  }
  else if(is.vector(data_t)){
    return(length(data_t))
  }
  else if(is.list(data_t)){
    return(lengths(data_t))
  }
  else if(is.atomic(data_t)){
    return(1)
  }
}

