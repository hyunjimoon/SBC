rhat_with_ignore <- function(dm, ignore_NAs = FALSE, ignore_inf = FALSE) {
  stopifnot(posterior::is_draws_matrix(dm))
  # if(ignore_NAs) {
  #   dm[is.na(dm)] <-
  # }
  posterior::rhat()
}
