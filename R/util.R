
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

if (link == 1) {
  p = brms:::inv_logit(eta)
} else if (link == 2) {
  p = dnorm(eta)
} else if (link == 3) {
  p = brms:::inv_cloglog(eta);
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
      param = brms:::inv_cloglog(eta)
    }
  param
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

intv_plot_save <- function(evolve_df, delivDir){
  for (v in names(evolve_df)){
    evolve_df_v <- evolve_df[[v]]
    intv <- subset_draws(mutate_variables(as_draws_df(lapply(evolve_df_v, as.numeric)), low1sd = (median + mad), up1sd = (median - mad)) , c("low1sd", "up1sd"))
    intv$iter <- as.numeric(rownames(evolve_df_v))
    intv <- reshape2::melt(intv, id.vars = "iter")
    intv <- filter(intv, variable == "low1sd" | variable == "up1sd")
    ggplot(intv, aes(x = as.numeric(iter), y = value,  color = variable) ) +
      geom_line() #+ ggtitle(sprintf("target par: %s, N: %s, M: %s ", pars, N, M))
    ggsave(file = file.path(delivDir, paste0(paste0(modelName, "_"), "evolove.png")), width = 5, height = 5)
  }
}

set_get_Dir <- function(modelName){
  scriptDir <- getwd()
  modDir <- file.path(scriptDir, "models")
  dir.create(file.path(scriptDir, "deliv"))
  delivDir <- file.path(scriptDir, "deliv", modelName)
  dir.create(delivDir)
  file <- file.path(modDir, paste0(modelName, ".stan"))
  mod <- cmdstan_model(file)
  return(list( mod = mod, modDir = modDir, file = file, delivDir = delivDir))
}

csv_save <- function(res, delivDir, type){
  if(type == "each") {
    write.csv(as_draws_df(res), file =  file.path(delivDir, paste0(paste0(cnt, "_"), "each.csv", sep = "")))
  } else if (type == "evolve"){
    for(v in names(res)) write.csv(res, file = file.path(paste0(paste0(delivDir, paste0(v, "_evolve_df.csv")))))
  } else if (type == "ecdf"){
    for(v in names(res)) write.csv(res, file =  file.path(delivDir, paste0(paste0(cnt, "_"), "ecdf.csv", sep = "")))
  } else if (type == "diagnositcs"){
    write.csv(res, file =  file.path(delivDir, "diagnositcs.csv"))
  }
}

pp_overlay_save <- function(param, param_next, cnt = 0, delivDir){
  g <- list()
  plotlist <- list()
  for (v in names(param)){
    g[[v]] <- bayesplot::ppc_dens_overlay(c(draws_of(param[[v]])), matrix(draws_of(param_next[[v]]), ncol = niterations(param[[v]])))
  }
  ggpubr::ggarrange(g[["a"]], nrow =length(names(param))) #TODO
  ggplot2::ggsave(file = file.path(delivDir, paste0(paste0(paste0(names(param), collapse = "", sep = "_"), cnt, "_"), "pp.png")),  bg = "white")
}
