#' A backend that wraps another backend with a caching scheme.
#'
#' The cache is stored in files named based on dataset and underlying
#' backend hashes, within the given directory, you can use
#' [cached_fit_filename()] to manually access the cache.
#'
#' The cached fit will try to re-emit all output, warnings and messages.
#' The ordering within each category is preserved, but the relative ordering of
#' e.g. outputs vs. messages (and all other combinations) may not be.
#'
#'
#' @export
#' @param cache_dir directory where the cache files (one file per fit) will
#' be stored
#' @param backend any other backend that does the actual computation
SBC_backend_cached <- function(cache_dir, backend) {
  if(!dir.exists(cache_dir)) {
    stop(paste0("cache_dir = '", cache_dir, "' is not an existing directory."))
  }
  structure(list(cache_dir = cache_dir,
                 backend = backend),
            class = "SBC_backend_cached")
}

#' Provide a file where a cached result of a fit will be/was stored.
#'
#' This is to allow you to directly access the cache, should you need to.
#' For normal fits, the cache file will always be a `list` readable by `readRDS` and have
#' element `$fit` for the actual fit. Other elements will contain the
#' text output, messages and warnings emitted by the fit.
#'
#' This is also used internally to cache some extra data (notably bridgesampling
#' results) using a non-empty `suffix` string.
#'
#' @export
cached_fit_filename <- function(cache_dir, backend, generated, suffix = "") {
  cache_basename <- paste0(
    "back_",
    SBC_backend_hash_for_cache(backend),
    "_ds_",
    rlang::hash(generated),
    suffix,
    ".rds"
  )
  cache_file <- file.path(cache_dir, cache_basename)
  return(cache_file)
}

#' @export
SBC_fit.SBC_backend_cached <- function(backend, generated, cores) {
  cache_file <- cached_fit_filename(backend$cache_dir, backend$backend, generated)
  if(file.exists(cache_file)) {
    message("Fit read from cache file '", cache_file, "'")
    res_and_output <- readRDS(cache_file)
    reemit_captured(res_and_output)
    fit <- SBC_backend_postprocess_cached_fit(backend$backend, generated, res_and_output$fit)
    return(fit)
  } else {
    res_and_output <- capture_all_outputs(
      SBC_fit(backend$backend, generated, cores)
    )
    names(res_and_output)[names(res_and_output) == "result"] <- "fit"

    reemit_captured(res_and_output)
    message("Storing fit in cache file '", cache_file, "'")
    res_and_output$fit <- SBC_backend_preprocess_fit_for_cache(backend$backend, generated, res_and_output$fit)
    saveRDS(res_and_output, file = cache_file)
    return(res_and_output$fit)
  }
}

#' @export
SBC_backend_hash_for_cache.SBC_backend_cached <- function(backend) {
  SBC_backend_hash_for_cache(backend$backend)
}

#' @export
SBC_backend_iid_draws.SBC_backend_cached <- function(backend) {
  SBC_backend_iid_draws(backend$backend)
}

#' @export
SBC_backend_default_thin_ranks.SBC_backend_cached <- function(backend) {
  SBC_backend_default_thin_ranks(backend$backend)
}

#' @export
SBC_fit_to_bridge_sampler.SBC_backend_cached <- function(backend, fit, generated, ...) {
  cache_file <- cached_fit_filename(backend$cache_dir, backend$backend, generated, suffix = "_bridgesampling")
  if(file.exists(cache_file)) {
    message("bridgesampler read from cache file '", cache_file, "'")
    res_and_output <- readRDS(cache_file)
    reemit_captured(res_and_output)
    return(res_and_output$bridgesampler)
  } else {
    res_and_output <- capture_all_outputs(
      SBC_fit_to_bridge_sampler(backend$backend, fit, generated, ...)
    )
    names(res_and_output)[names(res_and_output) == "result"] <- "bridgesampler"

    reemit_captured(res_and_output)
    message("Storing bridgesampler in cache file '", cache_file, "'")
    saveRDS(res_and_output, file = cache_file)
    return(res_and_output$bridgesampler)
  }
}

#' Allows the backend to do any pre-/post- processing on a fit stored to / loaded from cache.
#'
#' This is useful e.g. to provide a live `stanmodel` instance to an rstan fit,
#' as the one loaded from cache won't work for some cases.
#'
#' The deafult implementation does nothing.
#'
#' @export
SBC_backend_postprocess_cached_fit <- function(backend, generated, fit) {
  UseMethod("SBC_backend_postprocess_cached_fit")
}

#' @export
#' @rdname SBC_backend_postprocess_cached_fit
SBC_backend_preprocess_fit_for_cache <- function(backend, generated, fit) {
  UseMethod("SBC_backend_preprocess_fit_for_cache")
}


#' @export
SBC_backend_postprocess_cached_fit.default <- function(backend, generated, fit) {
  return(fit)
}

#' @export
SBC_backend_preprocess_fit_for_cache.default <- function(backend, generated, fit) {
  return(fit)
}
