

#' @description
#' Given a generator function, run the generator n times, and created draws_rvar for prior(parameters) and simulated data(generated)
#' Within each draws_rvar, n_chains is used as "n", which is the number of SBC iterations. For example, if we run SBC 10 times for 5 parameters,
#' the resulting draws_rvar will have dim = (1, 10, 5), with 10 chains and 5 parameters. For simulated data, if we run SBC 10 times and each
#' iteration yields 1000 posterior samples, the resulting draws_Rvar will have dim = (1000, 10, 1) with 1000 iterations and 10 chains.
#' Returned type is a named list of draws_rvars, for each parameters and generated
#'
#' @param generator The generator function
#' @param n_sbc_iterations the number of SBC iterations
#' @param ... Additional arguments to be passed to the generator function
#'
#' @return A named list with names = ("parameters", "generated"), which each contain draws_rvars corresponding to prior and simulated data
generator_to_draws_rvars <- function(generator, n_sbc_iterations, param, ...){
  merged_parameter_array <- NULL # a 3d array with dimensions(1, n_sbc_iterations, n_variables)
  merged_generated_array <- NULL # a 3d array with dim(n_generated_samples, n_sbc_iterations, 1)
  # we extend the array along the 2nd dimension, which is along chains

  # If there's an easier way to merge rvars along chains, then the double conversion wouldn't be needed
    for(iter in 1:n_sbc_iterations){
      generator_output <- do.call(generator, list(iter, ...))
      parameter_array <- posterior::as_draws_array(posterior::as_draws_rvars(generator_output[["parameters"]]))  # No clean way to directly transform generator output to draws_array
      generated_array <- posterior::as_draws_array(posterior::as_draws_rvars(list(y=posterior::rvar(generator_output[["generated"]]))))
      dimnames(parameter_array)[[2]] <- iter
      dimnames(generated_array)[[2]] <- iter
      if(iter == 1){
        merged_parameter_array <- parameter_array
        merged_generated_array <- generated_array
      }
      else{
        merged_parameter_array <- abind::abind(merged_parameter_array, parameter_array, along=2)
        merged_generated_array <- abind::abind(merged_generated_array, generated_array, along=2)
    }
  }
  list(parameters=posterior::as_draws_rvars(merged_parameter_array), generated=posterior::as_draws_rvars(merged_generated_array))
}


#' @description
#' Given a cmdstan model, a default data list, and a posterior::draws_rvars of simulated y samples, fit the model
#' n_sbc_iterations times each with data from simulated_y_draws_rvars, and return the merged posterior as a posterior::draws_rvars
#' object.
#'
#' @param cmdstan_model A CmdStanModel object
#' @param n_sbc_iterations Integer representing the number of iterations. Should equate to # of chains in sumulated_y_draws_rvars
#' @param data_list A list containing default data for the model. Note that "y" will be replaced with samples from simulate_y_draws_rvars
#' @param simulated_y_draws_rvars A posterior::draws_rvars containing simulated data. Should be of dim(n, n_sbc_iterations, 1)
#' @param ... Additional arguments to be passed to cmdstan_model$sample
#'
#' @return A posterior::draws_rvars object of dimension(n_fit_samples, n_sbc_iterations, n_variables) containing posterior samples
cmdstan_fits_to_draws_rvars <- function(cmdstan_model, n_sbc_iterations, data_list, simulated_y_draws_rvars, ...){
  if(n_sbc_iterations != posterior::nchains(simulated_y_draws_rvars)){
    stop("Number of chains in simulated_y_draws_rvars does not match n_sbc_iterations")
  }
  merged_fits_array <-NULL # a 3d array with dimensions(n_posterior_samples, n_sbc_iterations, n_variables)

  for(iter in 1:n_sbc_iterations){
    data_list[["y"]] <- as.vector(posterior::draws_of(posterior::subset_draws(simulated_y_draws_rvars, chain=iter)$y))
    model_fit <-do.call(cmdstan_model$sample, list(data=data_list, ...))
    fit_array <- posterior::as_draws_array(posterior::as_draws_rvars(model_fit$draws()))
    dimnames(fit_array)[[2]] <- iter

    if(iter == 1){
      merged_fits_array <- fit_array
    }
    else{
      merged_fits_array <- abind::abind(merged_fits_array, fit_array, along=2)
    }
  }
  posterior::as_draws_rvars(merged_fits_array)
}

#' @description
#' Given a posterior::draws_rvars types for both prior and posterior samples(refer to SBCWorkflow for specifications), create a
#' posterior::draws_rvars object which contains rank statistics (posterior < prior) for all parameters.
#'
#' @param prior_draws_rvars A posterior::draws_rvars type containing prior samples. Refer to SBCWorkflow for exact specifications
#' @param posteror_draws_rvars A posterior::draws_rvars type containing posterior samples. Refer to SBCWorkflow for exact specifications
#'
#' @return a Posterior::darws_rvars type containing rank statistcs. Dimensions are (1, n_sbc_iterations, n_variables)
calculate_rank_draws_rvars <- function(prior_draws_rvars, posterior_draws_rvars, param){
  n_vars <- posterior::nvariables(prior)
  sbc_iters <- posterior::nchains(prior)
  merged_rank_array <- NULL
  for(n_iter in 1:sbc_iters){
    rank_list <- list()
    for(n_var in 1:n_vars){
      var <- posterior::variables(prior)[[n_var]]

      prior_slice <- as.array(posterior::draws_of(posterior::subset_draws(prior, variable=var, chain=n_iter)[[var]]))
      iter_comp <- apply(posterior::draws_of(posterior::subset_draws(post, variable = var, chain=n_iter)[[var]]), 1, function(x){x < prior_slice})
      if(is.vector(iter_comp)){
        rank_list[[var]] <- sum(iter_comp, na.rm=TRUE)
      }
      else{
        rank_slice <- apply(iter_comp, c(1:(max(1,length(dim(iter_comp))-1))), function(x){sum(x, na.rm=TRUE)})
        # Need to test for scalar and matrices
        rank_list[[var]] <- rank_slice
      }

    }
    rank_array <- posterior::as_draws_array(posterior::as_draws_rvars(rank_list))
    dimnames(rank_array)[[2]] <- n_iter
    if(n_iter == 1){
      merged_rank_array <- rank_array
    }
    else{
      merged_rank_array <- abind::abind(merged_rank_array, rank_array, along=2)
    }
  }
  posterior::as_draws_rvars(merged_rank_array)
}



#############################################
#Legacy functions

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
            #mean_ <- summary_data[generate_index_name(parname, bitwOr(replace(tmp_indexes, accumulator, tmp_indexes[accumulator]+1), final_indexes)), "mean"]
            mean_ <- summary_data[generate_index_name(parname, bitwOr(replace(tmp_indexes, accumulator, tmp_indexes[accumulator]), final_indexes)), "mean"]
            stopifnot(!is.na(mean_))  # handle data.frames, which do not raise index errors, but returns NA
          },
          error = function(err){  # subscript out of bounds, backtrack and end component
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
}


#' @description
#' Given a named list, where indexes are parameter names with bracketed indexes(in the form of mu\[1\], mu\[2\]),
#' convert indexed parameters into actual R vectors and matrices.
#'
#' @param named_list A named list containing bracketed indexed parameter names as indexes.
#' @param parameter_dims a name list containing dimension info. Same as the output of infer_sequential_params(return_dim_info=TRUE)$dims
#'
#' @return A list with bracketed names removed, and with adequate vector/matrix structures.
par_list_to_structure = function(named_list, parameter_dims){
  unique_par_names <- unlist(unique(lapply(names(named_list), function(n){if(grepl("\\[", n)) sub("\\[.*", "", n) else n})))  # extract base parameter names(A[1, 2] -> A)
  generate_index_name <- function(name, ...){
    return(paste0(name, "[", paste0(c(...), collapse=","), "]"))  # receive indexes and return stan param name
  }
  return_list = list()
  for(i in 1:length(unique_par_names)){
    base_parname <- unique_par_names[[i]]

    dim_length <- length(parameter_dims[[base_parname]])
    par_dims <- parameter_dims[[base_parname]]


    if(dim_length == 1){
      if(par_dims == 1){  # scalar parameter
        return_list[[base_parname]] <-named_list[[base_parname]]
      }
      else{  # vector or 1d array
        return_list[[base_parname]] <- rep(0, parameter_dims[[i]])
        for(comp_index in 1:(parameter_dims[[i]])){
          return_list[[base_parname]][[comp_index]] <- named_list[[generate_index_name(base_parname, c(comp_index))]]
        }
      }
    }
    else if(dim_length > 1){
      n_components <- prod(par_dims)
      linear_array <- rep(NA, n_components)
      for(i in 1:n_components){
        arr_ind <- arrayInd(i, par_dims)[1, ]
        linear_array[[i]] <- as.numeric(named_list[[generate_index_name(base_parname, arr_ind)]])
      }
      return_list[[base_parname]] <- array(as.numeric(unlist(linear_array)), dim=par_dims)
    }
  }
  return(return_list)
}

#' @description
#' Given a named list, where elements are r data types, convert vectors and matrices to indexed element names.
#' For example, a vector \code{x} with 3 elements would be transformed into \code{x\[1\]}, \code{x\[2\]}, \code{x\[3\]}
#' Currently supports scalar, vector, and array types
#'
#' @param named_list A named list containing r elements
#'
#' @return A named list with indexed element values
decompose_structure_to_par_list <- function(named_list){
  generate_index_name <- function(name, ...){
    return(paste0(name, "[", paste0(c(...), collapse=","), "]"))
  }
  list_length <- length(named_list)
  list_names <- names(named_list)
  return_list <- list()
  for(n in 1:list_length){
    current_val <- named_list[[list_names[[n]]]]
    if(is.array(current_val)){
      linear_array <- as.vector(current_val)
      for(m in 1:prod(dim(current_val))){
        return_list[paste0(list_names[[n]], "[", paste(arrayInd(m, dim(current_val)), collapse=","), "]")] <- linear_array[m]
      }
    }
    else if(length(current_val) > 1){
      for(m in 1:length(current_val)){
        return_list[paste0(list_names[[n]], "[", m, "]")] <- current_val[m]
      }
    }
    else{
      return_list[list_names[[n]]] <- current_val
    }
  }
  return(return_list)
}

access_element_by_index <- function(data, index_list){
  if(length(index_list) == 0){
    return(data)
  }
  else if(length(index_list) == 1){
    return(data[[index_list[[1]]]])
  }
  else{
    return(cbind(unlist(index_list)))
  }
}
