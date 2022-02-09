#' Construct a generator used in the examples.
#'
#' @param example name of example
#' @param N size of the dataset the generator should simulate
#' @return an object that can be passed to [generate_datasets()]
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

#' Print the Stan code of a model used in the examples.
#'
#' @param example name of the example model.
#' @export
SBC_print_example_model <- function(example = c("normal_sd", "normal_var")) {
  example <- match.arg(example)

  example_program <- paste0(example, ".stan")
  code <- readLines(system.file(example_program, package = "SBC"))
  cat(code, sep = "\n")
}


#' Construct a backend to be used in the examples.
#'
#' @param example name of the example model
#' @param interface name of the interface to be used to fit the model
#' @export
SBC_example_backend <- function(example = c("normal_sd", "normal_var"),
                                interface = c("rstan", "cmdstanr")) {

  example <- match.arg(example)
  interface <- match.arg(interface)

  example_program <- paste0(example, ".stan")

  tmp <- file.path(tempdir(), example_program)
  if (!file.exists(tmp)) {
    file.copy(system.file(example_program, package = "SBC"), tmp)
  }

  if(interface == "cmdstanr") {
    mod <- cmdstanr::cmdstan_model(tmp)
    SBC_backend_cmdstan_sample(mod, chains = 2, iter_warmup = 400)
  } else if(interface == "rstan") {
    mod <- rstan::stan_model(tmp)
    SBC_backend_rstan_sample(mod, chains = 2, iter = 1400, warmup = 400)
  }
}

#' Combine an example backend with an example generator to provide full
#' results that can be used to test other functions in the package.
#'
#' @param example - name of the example. `normal_ok` is an example
#' where the generator matches the model
#' (using the `normal` generator and `normal_sd` backend), while
#' `normal_bad` is an example with a mismatch between the generator and backend
#' that manifests in SBC (`normal_bad` combines the `normal` generator with
#'  `normal_var` backend).
#' @param interface name of the interface to be used for the backend
#' @param N number of datapoints to simulate from the generator for each simulation
#' @param n_sims number of simulations to perform
#' @export
SBC_example_results <- function(example = c("normal_ok", "normal_bad"),
                                interface = c("rstan", "cmdstanr"),
                                N = 100, n_sims = 50) {
  example <- match.arg(example)
  interface <- match.arg(interface)
  if(example == "normal_ok") {
    generator <- SBC_example_generator(example = "normal", N = N)
    backend <- SBC_example_backend(example = "normal_sd", interface = interface)
  } else if (example == "normal_bad") {
    generator <- SBC_example_generator(example = "normal", N = N)
    backend <- SBC_example_backend(example = "normal_var", interface = interface)
  } else {
    stop("Invalid example")
  }



  compute_SBC(
    generate_datasets(generator,n_sims = n_sims),
    backend
  )
}
