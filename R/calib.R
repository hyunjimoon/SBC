library(cmdstanr)
library(dplyr)
library(tidyverse)
library(reshape)

iter_gen_inf <- function(sbc_obj, prior_theta, par_names, N, M, data){ #fixed params as prior_phi and add clampedStan
  sampled_y = sbc_obj$sample_y_tilde(prior_theta, data=data)
  post_theta = sbc_obj$sample_theta_bar_y(sampled_y, data=data, pars=par_names, fit_iter=M)
  NMP_tb <- as_tibble(post_theta)
  NMP_G <- NMP_tb %>%  gather(variable)
  NMP_G$par <- gsub("\\..*","",NMP_G$variable)
  NMP_G$prior <- sapply(NMP_G$variable, function(x){prior_theta[as.integer(gsub(".*\\.","",x)), gsub("\\..*","",x)]})
  Epost_bar_pri <- NMP_G %>% group_by(par, prior) %>% summarise('mean'  = round(mean(value),3))
  Varpost_bar_pri<- NMP_G %>% group_by(par, prior) %>% summarise('sd'  = sd(value))
  postVar <- inner_join(Epost_bar_pri %>% group_by(par) %>% summarise('VE' = sd(mean)), Varpost_bar_pri %>% group_by(par) %>% summarise('EV' = mean(sd^2))) %>%
    mutate("postVar" = VE + EV) %>% select(par, postVar)
  priVar = NMP_G %>% select(prior, par)  %>% group_by(par) %>% summarise('priVar' = sd(prior)^2)
  ratioVar = inner_join(postVar, priVar) %>% mutate("ratioVar" = postVar/priVar)
  return(ratioVar)
  }
