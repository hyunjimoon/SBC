
#' Combine two named lists and overwrite elements with the same name
#' using the value from args2
combine_args <- function(args1, args2) {
  if(is.null(names(args1)) || is.null(names(args2))) {
    c(args1, args2)
  } else {
    shared <- intersect(names(args1), names(args2))
    shared <- setdiff(shared, "")
    for(s in shared) {
      args1[[s]] <- args2[[s]]
    }
    c(args1, args2[!(names(args2) %in% shared)])
  }
}

##' Transform parameters from constrained to uncontrained
##'
##' @param param `draws_rvars` type parameter values
##' @return list of uncontrained parameters and transformation type
##' @export
tf_param <- function(param){
  tf <- list()
  for (tv in names(param)){
    if(all(param[[tv]] > 0) & all(param[[tv]] < 1)){
      tf[[tv]] <- "logit"
      param[[tv]] <- gtools::logit(param[[tv]])
    }else if(all(param[[tv]] > 0)){
      tf[[tv]] <- "log"
      param[[tv]] <- log(param[[tv]])
    }
  }
  return (list(param = param, tf = tf))
}

##' Inverse transform parameters from uncontrained to constrained
##'
##' @param param uncontrained parameter `draws_rvars`
##' @param tf list of transformation type
##' @return contrained parameter `draws_rvars`
##' @export
invtf_param <- function(param, tf){
  for (tv in names(param)){
    if(is.null(tf[[tv]])){
      param[[tv]] <- param[[tv]]
    }
    else if(tf[[tv]] == "logit"){
      param[[tv]] <- gtools::inv.logit(param[[tv]])
    }else if(tf[[tv]] == "log"){
      param[[tv]] <- exp(param[[tv]])
    }
    return (param)
  }
}


##' Transform parameters from constrained to uncontrained
##'
##' @param param a vector
##' @param tf string indicating transformation type
##' @return list containing uncontrained parameters and transformation type
##' @export
tf_param_vec <- function(param, tf){
  if(is.null(tf) || missing(tf)){
    param <- param
  }
  else if(tf == "logit"){
    param <- gtools::logit(param)
  }else if(tf == "log"){
    param <- log(param)
  }
  return (param)
}

##' Inverse transform parameters from uncontrained to constrained
##'
##' @param param a vector
##' @param link_type int indicating link type
##' @return constrained parameter vector
##' @export
invtf_param_vec <- function(param, link_type){
  if(is.null(link_type) || missing(link_type)){
    param <- param
  }
  else if(link_type == 1){
    param <- brms:::inv_logit(param)
  } else if (link_type == 2) {
    param = dnorm(param)
  } else if (link_type == 3) {
    param = brms:::inv_cloglog(eta);
    return (param)
  }
}


#'Maximal coupling of two univariate Normal distributions
#'from https://github.com/pierrejacob/debiasedhmc/blob/1a2eeeb041eea4e5c050e5188e7100f31e61e35b/R/gaussian_couplings.R
#'@description Sample from maximal coupling of two univariate Normal distributions,
#'specified through their means and standard deviations.
#'@param mu1 mean of first distribution
#'@param mu2 mean of second distribution
#'@param sigma1 standard deviation of first distribution
#'@param sigma2 standard deviation of second distribution
#'
#'@export
rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
  x <- rnorm(1, mu1, sigma1)
  if (dnorm(x, mu1, sigma1, log = TRUE) + log(runif(1)) < dnorm(x, mu2, sigma2, log = TRUE)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rnorm(1, mu2, sigma2)
      reject <- (dnorm(y, mu2, sigma2, log = TRUE) + log(runif(1)) < dnorm(y, mu1, sigma1, log = TRUE))
    }
    return(c(x,y))
  }
}


SBC_error <- function(subclass, message, call = sys.call(-1), ...) {
  structure(
    class = c(subclass, "SBC_error", "error", "condition"),
    list(message = message, call = call),
    ...
  )
}


require_package_version <- function(package, version, purpose) {
  if(!requireNamespace(package, quietly = TRUE)) {
    stop(paste0("Using ", purpose, " requires the '", package, "' package"))
  }
  # Cannot use `versionCheck` of `requireNamespace` as that doesn't work when
  # the package is already loaded. Note that `packageVersion` and `package_version`
  # are completely different methods
  if(packageVersion(package) < package_version(version)) {
    stop(paste0("SBC requires ", package, " version >= ", version, ", please update to use SBC."))
  }
}

require_brms_version <- function(purpose) {
  require_package_version("brms", "2.16.1", purpose)
}

require_cmdstanr_version <- function(purpose) {
  require_package_version("cmdstanr", "0.4.0", purpose)
}
