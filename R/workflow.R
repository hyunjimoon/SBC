
CMDSTAN_MODEL_CLASS_NAME <- "CmdStanModel"
BRMSFIT_MODEL_CLASS_NAME <- "brmsfit"

#' R6 class representing a fluid SBC workflow, boasting a highly customizable pipeline and supports intermediate updates
#' @export
SBCWorkflow <- R6::R6Class("SBCWorkflow",
 public = list(
   #' @field cmdstan_model A CmdStanModel object used for SBC
   #' @field brms_model_fit A brmsfit object used for SBC
   #' @field sim_function If using CmdStanModel, a simulator function, which should return simulated prior and data samples for SBC. Please refer \code{SBCWorkflow$initialize()} for details.
   #' @field calculated_ranks A named list of vectors which contained calculated ranks, or NULL if not calculated.
   #' @field prior_samples A posterior::draws_rvars of dimension(n_iterations=1, n_chains=n_sbc_iterations, n_variables=n_variables) which stores prior samples
   #' @field simulated_y A posterior::draws_rvars of dimension(n_iterations=n_simulated_y, n_chains=n_sbc_iterations, n_variables=1), which stores simulated y as "y"
   #' @field posterior_samples A posterior::draws_Rvars of dimension(n_iterations=n_posterior_samples, n_chains=n_sbc_iterations, n_variables=n_variables), which stores fitted posterior samples
   #' @field thin_factor Non-negative integer representing posterior thinning interval
   cmdstan_model = NULL,
   brms_model_fit = NULL,
   sim_function = NULL,
   calculated_ranks = NULL,
   prior_samples = NULL,  #  written on simulate
   simulated_y = NULL,  #  written on simulate
   posterior_samples = NULL, # written on fit_model
   thin_factor = NULL, # type is integer. written on initialize and fit_model

   #' R6 Initializer for SBCWorkflow
   #'
   #' @param model_obj A CmdStanModel or brmsfit object to use for fitting data.
   #' @param sim_function If using CmdStanModel, a simulator function. The function must return a
   #'   named list with elements \code{parameters} and \code{generated}.
   #'   Parameters should be a named list containing prior samplers for
   #'   parameters. Generated must be a 1 dimensional vector containing
   #'   simulated y data generated from the \code{parameters} samples.
   #'
   #' @return None
   initialize = function(model_obj, sim_function){
     if(!is.na(match(CMDSTAN_MODEL_CLASS_NAME, class(model_obj)))){
       if(missing(sim_function)){
         stop("If a cmdstan model is given, sim_function must also be supplied.")
       }
       self$cmdstan_model <- model_obj
       self$sim_function <- sim_function
     }
     else if(!is.na(match(BRMSFIT_MODEL_CLASS_NAME, class(model_obj)))){
       if(!missing(sim_function)){
         warning("A brm model has been supplied, but sim_function was also included in arguments. It will be ignored.")
       }
       self$brms_model_fit <- model_obj
     }
     else{
       stop(paste("The specified model object is neither a", CMDSTAN_MODEL_CLASS_NAME, "or a", BRMSFIT_MODEL_CLASS_NAME, ". Please recheck your model object."))
     }
   },

   #' Sample \eqn{\tilde{\theta} \sim P(\theta)} and \eqn{\tilde{y} \sim P(\tilde{y}) | \tilde{\theta} using \code{sim_function}
   #'
   #' @param n_sbc_iterations Number of prior sample sets to generate. Equates to number of SBC iterations
   #' @param ... Additional arguments to be passed to \code{sim_function}
   simulate = function(n_sbc_iterations, param = TRUE,  ...){  # number of simulation draws should be at least 1000
     if(!is.null(self$brms_model_fit)){
       #prior_sampler <- update(self$model_obj, sample_prior="only", chains=1, iter=n * 2, warmup = n)
     }
     else if(!is.null(self$cmdstan_model)){
       draws_rvars_list <- do.call(generator_to_draws_rvars, list(self$sim_function, n_sbc_iterations, param = param, ...))
       self$prior_samples <- draws_rvars_list[["parameters"]]
       self$simulated_y <- draws_rvars_list[["generated"]]
     }
     else{
       stop("No valid model specified for the workflow. Please check that your model is a valid brmsfit or a CmdStanModel object.")
     }
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
     # if(posterior::nchains(self$prior_samples) != posterior::nchains(self$simulated_y)){
     #   stop("The number of simulated prior and data do not match. Please rerun SBCWorkflow$simulate() and if using CmdStan, make sure your generator is well formed.")
     # }


     if(!is.null(self$brms_model_fit)){
       #brms
     }
     else if(!is.null(self$cmdstan_model)){
       sbc_iterations <- posterior::nchains(self$prior_samples)
       posterior_draws_rvars <- cmdstan_fits_to_draws_rvars(self$cmdstan_model, sbc_iterations, data, self$simulated_y, iter_sampling=sample_iterations * thin_factor, iter_warmup = warmup_iterations, chains=1, thin=thin_factor)
       self$posterior_samples <- posterior_draws_rvars
     }

     posterior_draws_rvars
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
       param <- posterior::variables(self$prior_samples)
     }
     self$calculated_ranks <- calculate_rank_rvars(self$prior_samples, self$posterior_samples, param = param)
     self$calculated_ranks

     # ranks <- array(rep(0, n_iters * n_params), dim=c(n_iters, n_params))
     # dimnames(ranks)[2] <- list(names(decomposed_prior))
     #
     # for(i in 1:n_iters){
     #   prior <- self$prior_samples[[i]]
     #   decomposed_prior <- decompose_structure_to_par_list(prior)
     #
     #   for(j in 1:n_params){
     #     regex_res <- gregexpr("(?<=\\[)(.+)(?=\\])", param[[j]], perl = TRUE)[[1]]
     #     capture_start <- attr(regex_res, "capture.start")
     #     capture_length <- attr(regex_res, "capture.length")
     #     if(capture_start != -1){
     #       index <- as.integer(unlist(substr(param[[j]], capture_start, capture_start + capture_length - 1)))
     #       ranks[i, param[[j]]] <- sum(self$posterior_samples[[i]][[param[[j]]]] < access_element_by_index(self$prior_samples[[i]][[base_param_names[[j]]]], index))
     #     }
     #     else{
     #       t1 <- self$posterior_samples[[i]][[param[[j]]]]
     #       t2 <- self$prior_samples[[i]][[param[[j]]]]
     #       #ranks[i, param[[j]]] <- sum(self$posterior_samples[[i]][[param[[j]]]] < self$prior_samples[[i]][[param[[j]]]])
     #       ranks[i, param[[j]]] <- sum(t1 < t2)
     #     }
     #   }
     # }
     # return(ranks)
   },
   add_steps = function(N){

   }
 )
)
