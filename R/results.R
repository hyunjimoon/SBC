#' @export
SBC_results <- function(ranks, z_scores, sds, fits) {
  #TODO argument validation
  structure(list(ranks = ranks, z_scores = z_scores, sds = sds, fits = fits), class = "SBC_results")
}

#' @export
compute_results <- function(datasets, backend, cores = getOption("mc.cores", 1), thin_ranks = 1) {
  #TODO use future for multiprocessing
  #TODO allow discarding fits after summarising to save memory
  fits <- list()
  warned_vars <- FALSE
  for(i in 1:length(datasets)) {
    fits[[i]] <- SBC_fit(backend, datasets$generated[[i]], cores = cores)
    fit_matrix <- SBC_fit_to_draws_matrix(fits[[i]])
    missing_vars <- setdiff(posterior::variables(datasets$parameters), posterior::variables(fit_matrix))
    if(length(missing_vars) > 0 && !warned_vars) {
      warning("Some variables missing in the fit: ", missing_vars)
    }

    shared_vars <- intersect(posterior::variables(datasets$parameters),
                             posterior::variables(fit_matrix))

    if(i == 1) {
      z_scores <- sds <- matrix(NA_real_, nrow = length(datasets), ncol = length(shared_vars))
      ranks <- matrix(NA_integer_, nrow = length(datasets), ncol = length(shared_vars))
    } else if(length(shared_vars) != ncol(ranks)) {
      stop("Differing number of variables across fits")
    }


    fit_matrix <- posterior::subset_draws(fit_matrix, variable = shared_vars)

    fit_thinned <- posterior::thin_draws(fit_matrix, thin_ranks)

    params <- posterior::subset_draws(datasets$parameters,
                                      draw = i,
                                      variable = shared_vars)

    last_ranks <- calculate_ranks_draws_matrix(params, fit_thinned)
    if(i == 1) {
      max_rank <- attr(last_ranks, "max_rank")
    } else {
      if(max_rank != attr(last_ranks, "max_rank")) {
        stop("Differening max_rank across fits")
      }
    }

    ranks[i, ] <- last_ranks
    #TODO: z-scores, sds
  }

  attr(ranks, "max_rank") <- max_rank
  colnames(ranks) <- shared_vars

  SBC_results(ranks = ranks, z_scores = z_scores, sds = sds, fits = fits)
}

#' Discrete uniform distribution allowing for varying lower and upper bounds.
#'
#' Internal, should not be exported.
#' Based on https://stats.stackexchange.com/a/3939/73129
rdunif <- function(n, a, b) {
  ceiling(runif(n, min = a - 1, max= b))
}

#' @param dm draws_matrix of the fit (assumed to be already thinned if that was necessary)
#' @export
calculate_ranks_draws_matrix <- function(params, dm) {
  #TODO validate input
  max_rank <- posterior::ndraws(dm)

  less_matrix <- sweep(dm, MARGIN = 2, STATS = params, FUN = "<")
  rank_min <- colSums(less_matrix)

  # When there are ties (e.g. for discrete parameters), the rank is currently drawn stochastically
  # among the ties
  equal_matrix <- sweep(dm, MARGIN = 2, STATS = params, FUN = "==")
  rank_range <- colSums(equal_matrix)

  ranks <- rank_min + rdunif(posterior::nvariables(dm), a = 0, b = rank_range)

  attr(ranks, "max_rank") <- max_rank
  ranks
}
