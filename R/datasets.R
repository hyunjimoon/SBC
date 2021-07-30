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
                   posterior::subset_draws(x$generated, draw = indices))
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
