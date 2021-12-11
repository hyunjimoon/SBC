
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
