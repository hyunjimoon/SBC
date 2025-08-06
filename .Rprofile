source("renv/activate.R")
source("~/.Rprofile")
# Allows to change how all vignettes are run at once (especially to test rstan)
options("SBC.vignettes_cmdstanr" = TRUE)

if (interactive() && requireNamespace("progressr", quietly = TRUE)) {
  options(progressr.handlers = progressr::handler_cli)

  ## Enable global progression updates
  if (getRversion() >= "4.0.0") progressr::handlers(global = TRUE)

}
