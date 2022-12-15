.onLoad <- function(libname, pkgname) {
  adjust_gamma <<- memoise::memoise(adjust_gamma, cache = cachem::cache_mem(max_size = 1024 * 32))
}
