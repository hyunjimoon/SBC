
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

list_of_values_to_draws_rvars <- function(list_of_values) {
    varnames <- names(list_of_values)
    if(is.null(varnames) || any(is.na(varnames)) ||
       any(varnames == "") || length(unique(varnames)) != length(varnames)) {
      stop(SBC_error("SBC_datasets_error", "All elements of list_of_values must have a unique name"))
    }

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
      purrr::map(list_of_values,
                 ~ posterior::rvar(array(.x, dim = c(1, guess_dims(.x))), dimnames = guess_dimnames(.x))
                 )
      )

    return(vars_rvars)
}
