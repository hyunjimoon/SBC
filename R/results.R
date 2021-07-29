SBC_results <- function(ranks, z_scores, sds, fits) {
  #TODO argument validation
  structure(list(ranks = ranks, z_scores = z_scores, sds = sds, fits = fits), class = "SBC_results")
}

compute_results <- function(datasets, backend, cores = getOption("mc.cores", 1), thin_ranks = 10) {
  #TODO use future for multiprocessing
  #TODO allow discarding fits after summarising to save memory
  fits <- list()
  z_scores <- sds <- matrix(NA_real_, nrow = length(datasets), ncol = nvariables(datasets$parameters))
  ranks <- matrix(NA_integer_, nrow = length(datasets), ncol = nvariables(datasets$parameters))
  warned_vars <- FALSE
  for(i in 1:length(datasets)) {
    fits[[i]] <- SBC_fit(backend, subset_draws(datasets$generated, draw = i))
    # TODO consider just using as_draws and not insist on a specific format until needed
    fit_rvars <- SBC_fit_to_draws_rvars(fits[[i]])
    missing_vars <- setdiff(variables(datasets$parameters), variables(fit_rvars))
    if(length(missing_vars) > 0 && !warned_vars) {
      warning("Some variables missing in the fit: ", missing_vars)
    }
    fit_thinned <- thin_draws(fit_rvars, thin_ranks)
    # TODO: continue here
  }

  SBC_results(ranks = ranks, z_scores = z_scores, sds = sds, fits = fits)
}
