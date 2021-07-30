new_SBC_datasets <- function(parameters, generated) {


  structure(list(parameters = parameters,
                 generated = generated),
            class = "SBC_datasets")
}

#' @export
validate_SBC_datasets <- function(x) {
  stopifnot(is.list(x))
  stopifnot(inherits(x, "SBC_datasets"))
  if(!posterior::is_draws_rvars(x$parameters)) {
    stop("SBC_datasets object has to have a 'parameters' field of type draws_rvars")
  }

  if(!is.list(x$generated)) {
    stop("SBC_datasets object has to have a 'generated' field of type list")
  }

  if(posterior::nchains(x$parameters) != 1) {
    stop("Needs one chain")
  }

  if(posterior::ndraws(x$parameters) != length(x$generated)) {
    stop("Needs equal no. of draws for parameters and length of generated")
  }

  x
}

#' Create new datasets object
#' @export
SBC_datasets <- function(parameters, generated) {
  x <-  new_SBC_datasets(parameters, generated)
  validate_SBC_datasets(x)
  x
}

#' @export
length.SBC_datasets <- function(x) {
  validate_SBC_datasets(x)
  posterior::ndraws(x$parameters)
}

#' @export
`[.SBC_datasets` <- function(x, indices) {
  validate_SBC_datasets(x)
  new_SBC_datasets(posterior::subset_draws(x$parameters, draw = indices),
                   x$generated[indices])
}


#' Generate datasets.
#'
#' @return object of class SBC_datasets
#' TODO: seed
#' @export
generate_datasets <- function(generator, n_datasets) {
  UseMethod("generate_datasets")
}

#'@export
list_function_SBC_generator <- function(f, ...) {
  stopifnot(is.function(f))
  structure(list(f = f, args = list(...)), class = "list_function_SBC_generator")
}

#'@export
generate_datasets.list_function_SBC_generator <- function(generator, n_datasets) {
  parameters_list <- list()
  generated <- list()
  for(iter in 1:n_datasets){
    generator_output <- do.call(generator$f, generator$args)
    #TODO check valid output
    parameters_list[[iter]] <- posterior::as_draws_rvars(generator_output$parameters)
    generated[[iter]] <- generator_output$generated
  }

  parameters <- do.call(posterior::bind_draws, args = c(parameters_list, list(along = "draw")))

  SBC_datasets(parameters, generated)
}

#'@export
function_SBC_generator <- function(f, ...) {
  stopifnot(is.function(f))
  structure(list(f = f, args = list(...)), class = "function_SBC_generator")
}

#'@export
generate_datasets.function_SBC_generator <- function(generator, n_datasets) {
  # TODO: check correct output
  do.call(generator$f, c(list(n_datasets = n_datasets), generator$args))
}

#'@export
brms_SBC_generator <- function(...) {
  model_data <- brms::make_standata(..., sample_prior = "only")
  class(model_data) <- NULL

  args <- list(...)

  if(!is.null(args$algorithm) && args$algorithm != "sampling") {
    stop("Algorithms other than sampling not supported yet")
  }

  compiled_model <- stanmodel_for_brms(...)



  structure(list(
    model_data = model_data,
    args = args,
    compiled_model = compiled_model),
    class = "brms_SBC_generator")
}

#' @export
generate_datasets.brms_SBC_generator <- function(generator, n_datasets) {
  #TODO pass args for chains, .... to sampling
  prior_fit <- generator$compiled_model$sample(data = generator$model_data,
                        iter_warmup = 1000, iter_sampling = n_datasets,
                        chains = 1)

  prior_fit_brms <- brmsfit_from_stanfit(prior_fit, generator$args)

  all_predictions <- brms::posterior_predict(prior_fit_brms)

  original_data <- generator$args$data
  processed_formula <- prior_fit_brms$formula

  generated <- list()
  log_likelihoods <- numeric(n_datasets)
  for(i in 1:n_datasets) {

    new_dataset <- original_data
    if(inherits(processed_formula, "brmsformula")) {
      new_dataset[[processed_formula$resp]] <- all_predictions[i, ]
    } else if (inherits(processed_formula, "mvbrmsformula")) {
      stop("Multivariate formulas not supported yet")
    } else {
      stop("Unrecognized formula")
    }
    generated[[i]] <- new_dataset

    ## Compute the likelihoods (observation model)
    ll <- log_lik(prior_fit_brms, newdata = new_dataset, subset = i)
    log_likelihoods[i] <- sum(ll)
  }

  # Add the observation likelihoods to the lp__ of the prior to
  # (hopefully) get total likelihood
  draws <- posterior::as_draws_rvars(prior_fit_brms$fit)
  draws$lp__ <- draws$lp__ + posterior::rvar(log_likelihoods, dim = 1)

  SBC_datasets(draws, generated)

}

draws_rvars_to_standata <- function(x) {
  stopifnot(posterior::ndraws(x) == 1)
  lapply(x, FUN = function(x_rvar) {
    res <- draws_of(x_rvar, with_chains = FALSE)
    #TODO figure out how to distinguish between scalar and array of size 1
    if(identical(dim(x_rvar), 1L)) {
      as.numeric(res[1,])
    } else {
      dim(res) <- dim(res)[2:length(dim(res))]
      res
    }
  })
}

