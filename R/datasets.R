new_SBC_datasets <- function(parameters, generated) {


  structure(list(parameters = parameters,
                 generated = generated),
            class = "SBC_datasets")
}

#' @export
validate_SBC_datasets <- function(x) {
  stopifnot(is.list(x))
  stopifnot(inherits(x, "SBC_datasets"))
  if(!posterior::is_draws_matrix(x$parameters)) {
    stop("SBC_datasets object has to have a 'parameters' field of type draws_matrix")
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

#' Combine multiple datasets together.
#' @param ... datasets to bind
#' @export
bind_datasets <- function(...) {
  args <- list(...)

  purrr::walk(args, validate_SBC_datasets)
  #TODO check identical par names

  parameters_list <- purrr::map(args, function(x) x$parameters)
  generated_list <- purrr::map(args, function(x) x$generated)

  new_SBC_datasets(do.call(posterior::bind_draws, parameters_list),
                   do.call(c, generated_list))
}

#' Generate datasets.
#'
#' @return object of class SBC_datasets
#' TODO: seed
#' @export
generate_datasets <- function(generator, n_datasets) {
  UseMethod("generate_datasets")
}

#' Generate datasets via a function that creates a single dataset.
#'
#' @param f function returning a list with elements `parameters`
#' (prior draws, a list or anything that can be converted to `draws_rvars`) and
#' `generated` (observed dataset, ready to be passed to backend)
#' @param ... Additional arguments passed to `f`
#'@export
function_SBC_generator <- function(f, ...) {
  stopifnot(is.function(f))
  structure(list(f = f, args = list(...)), class = "function_SBC_generator")
}


#' @export
generate_datasets.function_SBC_generator <- function(generator, n_datasets) {
  parameters_list <- list()
  generated <- list()
  for(iter in 1:n_datasets){
    generator_output <- do.call(generator$f, generator$args)
    if(!is.list(generator_output) ||
       is.null(generator_output$parameters) ||
       is.null(generator_output$generated)) {
      stop("The generating function has to return a list with elements `parameters`
      (that can be converted to `draws_rvars`) `generated`")
    }
    # TODO add a validate_input generic that would let backends impose additional checks
    # on generated data.

    # Directly converting to draws_matrix does not preserve arrays
    parameters_list[[iter]] <- posterior::as_draws_matrix(
      posterior::as_draws_rvars(generator_output$parameters))
    if(posterior::ndraws(parameters_list[[iter]]) != 1) {
      stop("The `parameters` element of the generated data must contain only
      a single draw")
    }
    generated[[iter]] <- generator_output$generated
  }

  parameters <- do.call(posterior::bind_draws, args = c(parameters_list, list(along = "draw")))

  SBC_datasets(parameters, generated)
}

#' Wrap a function the creates a complete dataset.
#'
#' This creates a very thin wrapper around the function and can store additional
#' arguments, but does not do anything more..
#'
#' Running:
#'
#' ```r
#' gen <- custom_SBC_generator(f, <<some other args>>)
#' datasets <- generate_datasets(gen, n_datasets = my_n_datasets)
#' ```
#'
#' is equivalent to just running
#'
#' ```r
#' datasets <- f(<<some other args>>, n_datasets = my_n_datasets)
#' ```
#'
#' So whenever you control the code calling `generate_datasets`,
#' it usually makes more sense to just create an `SBC_datasets`
#' object directly and avoid using `custom_SBC_generator` and `generate_datasets` at all.
#' `custom_SBC_generator` can however be useful, when a code you
#' do not control calls `generate_datasets` for you and the
#' built-in generators do not provide you with enough flexibility.
#'
#'
#' @param f function accepting at least an `n_datasets` argument and returning
#' and `SBC_datasets` object
#' @param ... Additional arguments passed to `f`
#' @export
custom_SBC_generator <- function(f, ...) {
  stopifnot(is.function(f))
  structure(list(f = f, args = list(...)), class = "custom_SBC_generator")
}

#'@export
generate_datasets.custom_SBC_generator <- function(generator, n_datasets) {
  res <- do.call(generator$f, combine_args(generator$args, list(n_datasets = n_datasets)))
  res <- validate_SBC_datasets(res)
  res
}

#' Create a brms generator.
#'
#' Brms generator uses a brms model with `sample_prior = "only"` to generate
#' new datasets.
#'
#' @param ... arguments passed to `brms::brm`
#' @param generate_lp whether to compute the overall log-likelihood of the model
#' as an additional parameter. This can be somewhat computationally expensive,
#' but improves sensitivity of the SBC process.
#' @export
brms_SBC_generator <- function(..., generate_lp = TRUE) {
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
    generate_lp = generate_lp,
    compiled_model = compiled_model),
    class = "brms_SBC_generator")
}

#' @export
generate_datasets.brms_SBC_generator <- function(generator, n_datasets) {
  #TODO pass args for control, warmup, .... to sampling
  if(inherits(generator$compiled_model, "CmdStanModel")) {
      args_for_fitting <- translate_rstan_args_to_cmdstan(generator$args, include_unrecognized = FALSE)
      args_for_fitting$data <- generator$model_data

      if(is.null(args_for_fitting$chains)) {
        args_for_fitting$chains <- 1
      }
      if(is.null(args_for_fitting$thin)) {
        args_for_fitting$thin <- 1
      }


      args_for_fitting$iter_sampling <- ceiling(n_datasets / args_for_fitting$chains) * args_for_fitting$thin

      args_for_fitting
      prior_fit <- do.call(generator$compiled_model$sample,
                           args_for_fitting)


      # Change once https://github.com/stan-dev/cmdstanr/issues/205
      # is resolved
      summ <- prior_fit$summary() # Can trigger warnings for treedepth/divergent ...
      max_rhat <- max(summ$rhat)
      if(max_rhat > 1.01) {
        message("Warning: Some rhats are > 1.01 indicating the prior was not explored well.\n",
                "The highest rhat is ", round(max_rhat, 2)," for ", summ$variable[which.max(summ$rhat)],
                "\nConsider adding warmup iterations (via 'warmup' argument).")
      }
      min_ess <- min(summ$ess_bulk)
      if(min_ess < n_datasets / 2) {
        message("Warning: Bulk effective sample size for some parameters is less than half the number of datasets.\n",
                "The lowest ESS_bulk/n_datasets is ", round(min_ess / n_datasets, 2)," for ", summ$variable[which.min(summ$ess_bulk)],
                "\nConsider increased thinning  (via 'thin' argument) .")
      }

  } else if (inherits(generator$compiled_model, "stanmodel")) {
    args_to_pass <- c("thin", "warmup", "control", "refresh")
    args_for_fitting <- c(
      list(object = generator$compiled_model,
           data = generator$model_data,
           chains = 1
      ),
      generator$args[intersect(args_to_pass, names(generator$args))]
    )

    if(is.null(args_for_fitting$warmup)) {
      args_for_fitting$warmup <- 1000
    }
    if(is.null(args_for_fitting$chains)) {
      args_for_fitting$chains <- 1
    }
    if(is.null(args_for_fitting$thin)) {
      args_for_fitting$thin <- 1
    }

    args_for_fitting$iter <- args_for_fitting$warmup + ceiling(n_datasets / args_for_fitting$chains) * args_for_fitting$thin

    prior_fit <- do.call(rstan::sampling, args_for_fitting)

  } else {
    stop("Invalid generator$compiled_model")
  }

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
    if(generator$generate_lp) {
      ll <- log_lik(prior_fit_brms, newdata = new_dataset, subset = i, cores = 1)
      log_likelihoods[i] <- sum(ll)
    }
  }

  # # Alternative code - compute LP in parallel. Seems not worth the overhead
  # # for the models I tested so far
  # if(generator$generate_lp) {
  #   log_likelihoods <- future.apply::future_mapply(
  #     FUN = function(new_dataset, i) {
  #       ll <- brms::log_lik(prior_fit_brms, newdata = new_dataset, subset = i, cores = 1)
  #       sum(ll)
  #       },
  #     generated, 1:n_datasets,
  #     future.chunk.size = default_chunk_size(n_datasets))
  # }


  # Add the observation likelihoods to the lp__ of the prior to
  # (hopefully) get total likelihood
  draws <- posterior::as_draws_matrix(prior_fit_brms$fit)
  if(generator$generate_lp) {
    draws[,"lp__"] <- draws[,"lp__"] + log_likelihoods
  } else {
    new_variables <- setdiff(posterior::variables(draws), "lp__")
    draws <- posterior::subset_draws(draws, variable = new_variables)
  }

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
