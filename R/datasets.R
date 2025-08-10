new_SBC_datasets <- function(variables, generated, var_attributes = NULL, parameters = NULL) {
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(missing(variables)) {
      variables <- parameters
    }
  }

  unvalidated <- structure(list(variables = variables,
                 generated = generated,
                 var_attributes = var_attributes),
            class = "SBC_datasets")

  return(validate_SBC_datasets(unvalidated))
}

#' @export
validate_SBC_datasets <- function(x) {
  stopifnot(is.list(x))
  stopifnot(inherits(x, "SBC_datasets"))
  if(!is.null(x$parameters)) {
    warning("Encountered old version of datasets using `parameters`, which is deprecated, will rename to `variables`.")
    if(is.null(x$variables)) {
      x$variables <- x$parameters
    }
  }
  if(!posterior::is_draws_matrix(x$variables)) {
    stop("SBC_datasets object has to have a 'variables' field of type draws_matrix")
  }

  if(!is.list(x$generated)) {
    stop("SBC_datasets object has to have a 'generated' field of type list")
  }

  if(posterior::nchains(x$variables) != 1) {
    stop("The `variables` draws_matrix needs exactly one chain.")
  }

  if(posterior::ndraws(x$variables) != length(x$generated)) {
    stop("Needs equal no. of draws for variables and length of generated")
  }

  x$var_attributes <- validate_var_attributes(x$var_attributes)

  if(!is.null(x$var_attributes)) {
    varnames_no_arrays <- unique(variable_names_to_var_attributes_names(posterior::variables(x$variables)))
    unrecognized_names <- setdiff(names(x$var_attributes), varnames_no_arrays)
    if(length(unrecognized_names) > 0) {
      message(paste0(unrecognized_names, collapse = ", "))
      stop("Some names of elements of `var_attributes` do not correspond to a name of a variable.\n",
           "Note that arrays are treated as a single variable and attributes must match the name.\n",
           "without brackets.")
    }
  }


  x
}

#' Create new `SBC_datasets` object.
#'
#' In most cases, you may want to use `generate_datasets` to build the object, but
#' for full control, you can also create datasets directly via this function.
#'
#' @param variables draws of "true" values of unobserved parameters or other derived variables.
#' An object of class `draws_matrix` (from the `posterior` package)
#' @param generated a list of objects that can be passed as data to the backend you plan to use.
#' (e.g. list of values for Stan-based backends, a data frame for `SBC_backend_brms`)
#' @param parameters DEPRECATED. Use variables instead.
#' @export
SBC_datasets <- function(variables, generated, var_attributes = NULL, parameters = NULL) {
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(missing(variables)) {
      variables <- parameters
    }
  }
  x <- new_SBC_datasets(variables, generated, var_attributes)
  return(x)
}

#' @export
length.SBC_datasets <- function(x) {
  validate_SBC_datasets(x)
  posterior::ndraws(x$variables)
}

#' Subset an `SBC_datasets` object.
#'
#' @export
`[.SBC_datasets` <- function(x, indices) {
  validate_SBC_datasets(x)
  if(is.logical(indices)) {
    stopifnot(length(indices) == length(x))
    indices <- which(indices)
  }
  res <- new_SBC_datasets(posterior::subset_draws(x$variables, draw = indices, unique = FALSE),
                   x$generated[indices],
                   var_attributes = x$var_attributes)
  return(res)
}

#' Combine multiple datasets together.
#' @param ... datasets to bind
#' @export
bind_datasets <- function(...) {
  args <- list(...)

  if(length(args) == 0) {
    stop("Must bind at least 1 dataset")
  }

  purrr::walk(args, validate_SBC_datasets)
  #TODO check identical par names (throws error implicitly in bind_draws)

  variables_list <- purrr::map(args, function(x) x$variables)
  generated_list <- purrr::map(args, function(x) x$generated)

  var_attributes_list <- purrr::map(args, function(x) x$var_attributes)
  purrr::walk(var_attributes_list, \(attr) {
    if(!identical(attr, var_attributes_list[[1]])) {
      stop("All var_attributes of datasets to be bound must be identical")
    }
  })

  new_SBC_datasets(do.call(posterior::bind_draws, c(variables_list, list(along = "draw"))),
                   do.call(c, generated_list),
                   var_attributes = var_attributes_list[[1]])
}

#' Generate datasets.
#'
#' @param generator a generator object - build e.g. via `SBC_generator_function` or
#'  `SBC_generator_brms`.
#' @param n_sims the number of simulated datasets to use
#' @param n_datasets DEPRECATED, use `n_sims` instead.
#' @return object of class `SBC_datasets`
#' @export
generate_datasets <- function(generator, n_sims, n_datasets = NULL) {
  UseMethod("generate_datasets")
}

#' Generate datasets via a function that creates a single dataset.
#'
#'
#' The generator supports progress reporting with the
#' [`progressr`](https://progressr.futureverse.org/) package.
#' Run `progressr::handlers(global = TRUE)` to enable progress reporting
#' globally, see the `progressr` docs for more options.
#'
#'
#' @param f function returning a list with elements `variables`
#' (prior draws, a list or anything that can be converted to `draws_rvars`) and
#' `generated` (observed dataset, ready to be passed to backend)
#' @param future.chunk.size if finite, datasets will be generated in parallel using
#' the `future` package with this chunk size (see [future.apply::future_lapply]). Defaults to the `SBC.generator_chunk_size` option, which
#' defaults to +Inf (no paralellism).
#' @param ... Additional arguments passed to `f`
#'@export
SBC_generator_function <- function(f, future.chunk.size = getOption("SBC.generator_chunk_size", +Inf), future.globals = TRUE ,...) {
  stopifnot(is.function(f))
  structure(list(f = f, future.chunk.size = future.chunk.size, future.globals = future.globals, args = list(...)), class = "SBC_generator_function")
}


#' @export
generate_datasets.SBC_generator_function <- function(generator, n_sims, n_datasets = NULL) {
  if(!is.null(n_datasets)) {
    warning("n_datasets argument is deprecated, use n_sims instead")
    if(missing(n_sims)) {
      n_sims <- n_datasets
    }
  }

  stopifnot(n_sims >= 1)

  warned_parameters <- FALSE

  if(requireNamespace("progressr", quietly = TRUE)) {
    progressor <- progressr::progressor(n_sims)
  } else {
    progressor <- NULL
  }


  generate_single <- function() {
    generator_output <- do.call(generator$f, generator$args)
    # Ensuring backwards compatibility
    if(!is.null(generator_output$parameters)) {
      if(!warned_parameters) {
        warning("Generator function returns a list with element `parameters`, which is deprecated. Return `variables` instead.")
        warned_parameters <<- TRUE
      }
      if(is.null(generator_output$variables)) {
        generator_output$variables <- generator_output$parameters
      }
    }

    if(!is.list(generator_output) ||
       is.null(generator_output$variables) ||
       is.null(generator_output$generated)) {
      stop(SBC_error("SBC_datasets_error",
                     "The generating function has to return a list with elements `variables`
      (that can be converted to `draws_rvars`) and `generated`.
      Optionally a `var_attributes` element can also be returned."))
    }

    varnames <- names(generator_output$variables)
    if(is.null(varnames) || any(is.na(varnames)) ||
       any(varnames == "") || length(unique(varnames)) != length(varnames)) {
      stop(SBC_error("SBC_datasets_error", "All elements of $variables must have a unique name"))
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

    guess_dimnames <- function(x) {
      if(!is.null(dimnames(x))) {
        dimnames(x)
      } else if(!is.null(names(x))) {
        list(names(x))
      } else {
        NULL
      }
    }

    vars_rvars <-
      do.call(
        posterior::draws_rvars,
        purrr::map(generator_output$variables,
                   ~ posterior::rvar(array(.x, dim = c(1, guess_dims(.x))), dimnames = guess_dimnames(.x))
        )
      )

    vars_dm <- posterior::as_draws_matrix(vars_rvars)
    if(posterior::ndraws(vars_dm) != 1) {
      stop("The `variables` element of the generator output must contain only
      a single draw")
    }
    generator_output$variables <- vars_dm

    if(!is.null(progressor)) {
      progressor()
    }

    return(generator_output)
  }

  if(is.finite(generator$future.chunk.size)) {
    all_outputs <- future.apply::future_replicate(n_sims, generate_single(),
                                                  future.chunk.size = generator$future.chunk.size,
                                                  future.globals = generator$future.globals,
                                                  simplify = FALSE)
  } else {
    all_outputs <- replicate(n_sims, generate_single(), simplify = FALSE)
  }

  variables_dm_list <- purrr::map(all_outputs, \(x) x$variables)
  generated <- purrr::map(all_outputs, \(x) x$generated)

  # Process var_attributes
  var_attributes <- all_outputs[[1]]$var_attributes
  if(n_sims >= 2) {
    for (iter in 2:n_sims) {
      if (is.null(var_attributes) != is.null(all_outputs[[iter]]$var_attributes)) {
        stop("Some simulations returned $var_attributes, while others did not.")
      }
      if (!is.null(var_attributes) &&
          !identical(var_attributes, all_outputs[[iter]]$var_attributes)) {
        stop(
          "Some simulations returned different $var_attributes than others.
             Attributes must be identical for all simulations."
        )
      }
    }
  }

  variables <- do.call(posterior::bind_draws, args = c(variables_dm_list, list(along = "draw")))

  SBC_datasets(variables, generated, var_attributes = var_attributes)
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
#' datasets <- generate_datasets(gen, n_sims = my_n_sims)
#' ```
#'
#' is equivalent to just running
#'
#' ```r
#' datasets <- f(<<some other args>>, n_sims = my_n_sims)
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
#' @param f function accepting at least an `n_sims` argument and returning
#' and `SBC_datasets` object
#' @param ... Additional arguments passed to `f`
#' @export
SBC_generator_custom <- function(f, ...) {
  stopifnot(is.function(f))
  structure(list(f = f, args = list(...)), class = "SBC_generator_custom")
}

#'@export
generate_datasets.SBC_generator_custom <- function(generator, n_sims, n_datasets = NULL) {
  if(!is.null(n_datasets)) {
    warning("n_datasets argument is deprecated, use n_sims instead")
    if(missing(n_sims)) {
      n_sims <- n_datasets
    }
  }
  res <- do.call(generator$f, combine_args(generator$args, list(n_sims = n_sims)))
  res <- validate_SBC_datasets(res)
  res
}


#' @export
draws_rvars_to_standata <- function(x) {
  res <- list()
  for(i in 1:posterior::ndraws(x)) {
    # TODO use direct indexing - subset_draws is unnecessarily slow
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
    warning("Cannot reliably estimate prior_sd with less than 50 simulations.\n",
            "Note that you can generate extra simulations that you don't actually fit and use those to estimate prior sd.")
  }
  if(length(datasets) < 2) {
    stop("Cannot estimate prior sd with less than 2 simulations")
  }

  sds_df <- posterior::summarise_draws(datasets$variables, sd)
  sds_vec <- sds_df$sd
  names(sds_vec) <- sds_df$variable

  sds_vec
}
