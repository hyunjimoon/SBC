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

  expect_identical(posterior::variables(res$parameters), c("mu", "sigma"))
  expect_identical(names(res$generated[[1]]), c("y1", "y2"))

  expect_equal(posterior::ndraws(res$parameters), 7)

})

test_that("subsetting datasets", {
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

  res_subs <- res[3:5]

  expect_identical(res_subs$parameters, posterior::subset_draws(res$parameters, draw = 3:5))
  expect_identical(res_subs$generated[[1]], res$generated[[3]])
  expect_identical(res_subs$generated[[2]], res$generated[[4]])
  expect_identical(res_subs$generated[[3]], res$generated[[5]])
  expect_equal(length(res_subs$generated), 3)

})
