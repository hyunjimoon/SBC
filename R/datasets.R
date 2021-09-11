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

#' Create new `SBC_datasets` object.
#'
#' In most cases, you may want to use `generate_datasets` to build the object, but
#' for full control, you can also create datasets directly via this function.
#'
#' @param parameters samples of "true" values of unobserved parameters.
#' An object of class `draws_matrix` (from the `posterior` package)
#' @param generated a list of objects that can be passed as data to the backend you plan to use.
#' (e.g. list of values for Stan-based backends, a data frame for `SBC_backend_brms`)
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

#' Subset an `SBC_datasets` object.
#'
#' @export
`[.SBC_datasets` <- function(x, indices) {
  validate_SBC_datasets(x)
  new_SBC_datasets(posterior::subset_draws(x$parameters, draw = indices, unique = FALSE),
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

  new_SBC_datasets(do.call(posterior::bind_draws, c(parameters_list, list(along = "draw"))),
                   do.call(c, generated_list))
}

#' Generate datasets.
#'
#' @param generator a generator object - build e.g. via `SBC_generator_function` or
#'  `SBC_generator_brms`.
#' @return object of class `SBC_datasets`
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
SBC_generator_function <- function(f, ...) {
  stopifnot(is.function(f))
  structure(list(f = f, args = list(...)), class = "SBC_generator_function")
}


#' @export
generate_datasets.SBC_generator_function <- function(generator, n_datasets) {
  parameters_list <- list()
  generated <- list()
  for(iter in 1:n_datasets){
    generator_output <- do.call(generator$f, generator$args)
    if(!is.list(generator_output) ||
       is.null(generator_output$parameters) ||
       is.null(generator_output$generated)) {
      stop(SBC_error("SBC_datasets_error",
      "The generating function has to return a list with elements `parameters`
      (that can be converted to `draws_rvars`) `generated`"))
    }

    parnames <- names(generator_output$parameters)
    if(is.null(parnames) || any(is.na(parnames)) ||
       any(parnames == "") || length(unique(parnames)) != length(parnames)) {
      stop(SBC_error("SBC_datasets_error", "All elements of $parameters must have a unique name"))
    }
    # TODO add a validate_input generic that would let backends impose additional checks
    # on generated data.

    # Directly converting to draws_matrix does not preserve arrays
    guess_dims <- function(x) {
      if(!is.null(dim(x))) {
        dim(x)
      } else {
        if(length(x) > 1) {
          length(x)
        } else {
          NULL
        }
      }
    }

    params_rvars <-
      do.call(
      posterior::draws_rvars,
      purrr::map(generator_output$parameters,
                 ~ posterior::rvar(array(.x, dim = c(1, guess_dims(.x))))
                 )
      )
    parameters_list[[iter]] <- posterior::as_draws_matrix(params_rvars)
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
#' gen <- SBC_generator_custom(f, <<some other args>>)
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
#' object directly and avoid using `SBC_generator_custom` and `generate_datasets` at all.
#' `SBC_generator_custom` can however be useful, when a code you
#' do not control calls `generate_datasets` for you and the
#' built-in generators do not provide you with enough flexibility.
#'
#'
#' @param f function accepting at least an `n_datasets` argument and returning
#' and `SBC_datasets` object
#' @param ... Additional arguments passed to `f`
#' @export
SBC_generator_custom <- function(f, ...) {
  stopifnot(is.function(f))
  structure(list(f = f, args = list(...)), class = "SBC_generator_custom")
}

#'@export
generate_datasets.SBC_generator_custom <- function(generator, n_datasets) {
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
SBC_generator_brms <- function(..., generate_lp = TRUE) {
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
    class = "SBC_generator_brms")
}

#' @export
generate_datasets.SBC_generator_brms <- function(generator, n_datasets) {
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
                "The highest rhat is ", round(max_rhat, 2)," for ", summ$parameter[which.max(summ$rhat)],
                "\nConsider adding warmup iterations (via 'warmup' argument).")
      }
      min_ess <- min(summ$ess_bulk)
      if(min_ess < n_datasets / 2) {
        message("Warning: Bulk effective sample size for some parameters is less than half the number of datasets.\n",
                "The lowest ESS_bulk/n_datasets is ", round(min_ess / n_datasets, 2)," for ", summ$parameter[which.min(summ$ess_bulk)],
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

#' @export
draws_rvars_to_standata <- function(x) {
  res <- list()
  for(i in 1:posterior::ndraws(x)) {
    res[[i]] <- draws_rvars_to_standata_single(posterior::subset_draws(x, draw = i))
  }
  res
}

#' @export
draws_rvars_to_standata_single <- function(x) {
  stopifnot(posterior::ndraws(x) == 1)
  lapply(x, FUN = function(x_rvar) {
    res <- posterior::draws_of(x_rvar, with_chains = FALSE)
    #TODO figure out how to distinguish between scalar and array of size 1
    if(identical(dim(x_rvar), 1L)) {
      as.numeric(res[1,])
    } else {
      dim(res) <- dim(res)[2:length(dim(res))]
      res
    }
  })
}

#' Calculate  prior standard deviation of a dataset
#'
#' @returns a named vector of prior SDs
#' @export
calculate_prior_sd <- function(datasets) {
  validate_SBC_datasets(datasets)
  # TODO this is a hack - there has to be a better diagnostic to get whether
  # our sd estimate is good (probably via MCSE?)
  if(length(datasets) < 50) {
    warning("Cannot reliably estimate prior_sd with less than 50 datasets.\n",
            "Note that you can generate extra datasets that you don't actually fit and use those to estimate prior sd.")
  }
  if(length(datasets) < 2) {
    stop("Cannot estimate prior sd with less than 2 datasets")
  }

  sds_df <- posterior::summarise_draws(datasets$parameters, sd)
  sds_vec <- sds_df$sd
  names(sds_vec) <- sds_df$variable

  sds_vec
}
