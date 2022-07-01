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
SBC_print_example_model <- function(example = c("normal_sd", "normal_bad"),
                                    interface = c("rstan", "cmdstanr", "rjags")) {
  #Backward compatibility
  if(identical(example, "normal_var")) {
    example <- "normal_bad"
  }

  example <- match.arg(example)
  interface <- match.arg(interface)

  if(interface %in% c("rstan", "cmdstanr")) {
    example_program <- paste0(example, ".stan")
  } else if(interface == "rjags") {
    example_program <- paste0(example, ".jags")
  }
  code <- readLines(system.file(example_program, package = "SBC"))
  cat(code, sep = "\n")
}


#' Construct a backend to be used in the examples.
#'
#' Note that this will involve compiling a Stan model and may take a while.
#'
#' @param example name of the example model. `normal_sd` is a simple model fitting
#'  a normal distribution parametrized as mean and standard deviation.
#'  `normal_bad` is a model that _tries_ to implement the `normal_sd` model,
#'  but assumes an incorrect parametrization of the normal distribution.
#'  For Stan-based backends, the model is written as if Stan parametrized
#'  normal distribution with precision (while Stan uses sd), for JAGS-based
#'  backends the model is written as if JAGS parametrized normal distribution
#'  with sd (while JAGS uses precision).
#' @param interface name of the interface to be used to fit the model
#' @export
SBC_example_backend <- function(example = c("normal_sd", "normal_bad"),
                                interface = c("rstan", "cmdstanr", "rjags")) {

  #Backward compatibility
  if(identical(example, "normal_var")) {
    example <- "normal_bad"
  }

  example <- match.arg(example)
  interface <- match.arg(interface)

  if(interface %in% c("cmdstanr", "rstan")) {
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
  } else if(interface == "rjags") {
    model_file <- system.file(paste0(example, ".jags"), package = "SBC")
    SBC_backend_rjags(file = model_file, n.iter = 5000, n.burnin = 5000,
                      thin = 10, n.chains = 2,
                      variable.names = c("mu", "sigma"))
  } else {
    stop("Invalid interface")
  }
}

#' Combine an example backend with an example generator to provide full
#' results that can be used to test other functions in the package.
#'
#' Except for `example = "visualizations"`, all examples will actually
#' compile and fit Stan models and thus may take a while to complete.
#'
#' @param example - name of the example. `normal_ok` is an example
#' where the generator matches the model
#' (using the `normal` generator and `normal_sd` backend), while
#' `normal_bad` is an example with a mismatch between the generator and backend
#' that manifests in SBC (`normal_bad` combines the `normal` generator with
#'  `normal_bad` backend). `visualizations` creates a purely artificial results
#'  that are meant to showcase the built-in plots (the `interface` parameter will
#'  be ignored).
#' @param interface name of the interface to be used for the backend
#' @param N number of datapoints to simulate from the generator for each simulation
#' @param n_sims number of simulations to perform
#' @export
SBC_example_results <- function(example = c("normal_ok", "normal_bad", "visualizations"),
                                interface = c("rstan", "cmdstanr", "rjags"),
                                N = 100, n_sims = 50) {
  example <- match.arg(example)
  interface <- match.arg(interface)
  if(example == "normal_ok") {
    generator <- SBC_example_generator(example = "normal", N = N)
    backend <- SBC_example_backend(example = "normal_sd", interface = interface)
  } else if (example == "normal_bad") {
    generator <- SBC_example_generator(example = "normal", N = N)
    backend <- SBC_example_backend(example = "normal_bad", interface = interface)
  } else if (example == "visualizations") {

    df_x <- seq(-4, 4, length.out = 400)
    prior_df <- tidyr::crossing(data.frame(x = df_x, density = dnorm(df_x), type = "Correct"),
                         variable = c("Exact match",
                         "Model too certain",
                         "Model too uncertain",
                         "Model underestimating",
                         "Model overestimating",
                         "Some extra-low estimates"))

    generator <- SBC_generator_function(function() {
      list(
        variables = list(
          "Exact match" = rnorm(1),
          "Model too certain" = rnorm(1),
          "Model too uncertain" = rnorm(1),
          "Model underestimating" = rnorm(1),
          "Model overestimating" = rnorm(1),
          "Some extra-low estimates" = rnorm(1)
        ),
        generated = list()
      )
    })

    posterior_df <- rbind(
      data.frame(variable = "Exact match", x = df_x, density = dnorm(df_x)),
      data.frame(variable = "Model too certain", x = df_x, density = dnorm(df_x, sd = 1/3)),
      data.frame(variable = "Model too uncertain", x = df_x, density =  dnorm(df_x, sd = 2)),
      data.frame(variable = "Model underestimating", x = df_x, density =  dnorm(df_x, mean = -1)),
      data.frame(variable = "Model overestimating", x = df_x, density =  dnorm(df_x, mean = 1)),
      data.frame(variable = "Some extra-low estimates", x = df_x,
                 density =  0.1 * dnorm(df_x, mean = -3, sd = 0.1) + 0.9 * dnorm(df_x))
    )
    posterior_df$type = "Observed"

    backend <- SBC_backend_mock_rng(
      "Exact match" = ~ rnorm(.),
      "Model too certain" = ~ rnorm(., sd = 1/3),
      "Model too uncertain" = ~ rnorm(., sd = 2),
      "Model underestimating" = ~ rnorm(., mean = -1),
      "Model overestimating" = ~ rnorm(., mean = 1),
      "Some extra-low estimates" = function(.) { if(runif(1) < 0.1) { rnorm(., mean = -3, sd = 0.1) } else { rnorm(.) }},
      n_draws = 100
    )

    res <- compute_SBC(
      generate_datasets(generator,n_sims = n_sims),
      backend
    )

    attr(res, "density_df") <- rbind(prior_df, posterior_df)

    return(res)
  } else {
    stop("Invalid example")
  }



  compute_SBC(
    generate_datasets(generator,n_sims = n_sims),
    backend
  )
}
