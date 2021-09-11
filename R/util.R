
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

pp_overlay_save <- function(param, param_next, cnt = 0){
  g <- list()
  plotlist <- list()
  for (v in names(param)){
    g[[v]] <- ppc_dens_overlay(c(draws_of(param[[v]])), matrix(draws_of(param_next[[v]]), ncol = niterations(param[[v]])))
  }
  ggarrange(g[["loc"]], g[["scale"]], nrow =length(names(param))) #TODO
  ggsave(file = file.path(delivDir, paste0(paste0(paste0(names(param), collapse = "", sep = "_"), cnt, "_"), "pp.png")),  bg = "white")
}
