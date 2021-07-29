test_that("Generating datasets via functions", {
  list_function <- function(N) {
    mu <- rnorm(1)
    sigma <- abs(rnorm(1))
    y1 <- rnorm(N, mu, sigma)
    y2 <- rnorm(2 * N, mu + 5, sigma)
    list(parameters = list(mu = mu, sigma = sigma),
         generated = list(y1 = y1, y2 = y2))
  }

  res <- generate_datasets(
    list_function_SBC_generator(list_function, N = 10),
    n_datasets = 7)

  expect_true(length(res) == 7)

  expect_identical(variables(res$parameters), c("mu", "sigma"))
  expect_identical(variables(res$generated), c("y1", "y2"))

  expect_equal(length(res$parameters$mu), 1)
  expect_equal(length(res$parameters$sigma), 1)
  expect_equal(length(res$generated$y1) , 10)
  expect_equal(length(res$generated$y2) , 20)

  # The same, but with a direct function
  direct_function <- function(n_datasets, N) {
    mu <- posterior::rvar_rng(rnorm, 1, ndraws = n_datasets)
    sigma <- abs(posterior::rvar_rng(rnorm, 1, ndraws = n_datasets))

    y1 <- rvar_rng(rnorm, N, mu, sigma)
    y2 <- rvar_rng(rnorm, 2 * N, mu + 5, sigma)

    SBC_datasets(parameters = posterior::draws_rvars(mu = mu, sigma = sigma),
                 generated = posterior::draws_rvars(y1 = y1, y2 = y2))
  }

  res <- generate_datasets(
    function_SBC_generator(direct_function, N = 10),
    n_datasets = 7)

  expect_true(length(res) == 7)

  expect_identical(variables(res$parameters), c("mu", "sigma"))
  expect_identical(variables(res$generated), c("y1", "y2"))

  expect_equal(length(res$parameters$mu), 1)
  expect_equal(length(res$parameters$sigma), 1)
  expect_equal(length(res$generated$y1) , 10)
  expect_equal(length(res$generated$y2) , 20)

})

