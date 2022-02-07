#' @export
SBC_example_generator <- function(example = c("normal"), N = 100) {
  example <- match.arg(example)
  if(example == "normal") {
    generator_func <- function(N) {
      mu <- rnorm(1, 0, 1);
      sigma <- abs(rnorm(1, 0, 1))
      y <- rnorm(N, mu, sigma)
      list(
        variables = list(
          mu = mu,
          sigma = sigma
        ),
        generated = list(
          N = N,
          y = y
        )
      )
    }
  } else {
    stop("Invalid dataset example")
  }

  SBC_generator_function(generator_func, N = N)
}

#' @export
SBC_print_example_model <- function(example = c("normal_sd", "normal_var")) {
  example_program <- paste0(example, ".stan")
  code <- readLines(system.file(example_program, package = "SBC"))
  cat(code, sep = "\n")
}


#' @export
SBC_example_backend <- function(example = c("normal_sd", "normal_var"),
                                interface = c("cmdstanr", "rstan")) {
  example_program <- paste0(example, ".stan")

  tmp <- file.path(tempdir(), example_program)
  if (!file.exists(tmp)) {
    file.copy(system.file(example_program, package = "SBC"), tmp)
  }

  if(interface == "cmdstanr") {
    mod <- cmdstanr::cmdstan_model(tmp)
    SBC_backend_cmdstan_sample(mod, chains = 2, iter_warmup = 200)
  } else if(interface == "rstan") {
    mod <- rstan::stan_model(tmp)
    SBC_backend_rstan_sample(mod, chains = 2, iter = 1200, warmup = 200)
  }
}


SBC_example_results <- function(example = c("normal_pass", "normal_fail"),
                                interface = c("cmdstanr", "rstan")) {

}
