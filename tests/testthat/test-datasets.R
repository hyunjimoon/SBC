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
    SBC_generator_function(list_function, N = 10),
    n_datasets = 7)

  expect_true(length(res) == 7)

  expect_identical(posterior::variables(res$parameters), c("mu", "sigma"))
  expect_identical(names(res$generated[[1]]), c("y1", "y2"))

  expect_equal(posterior::ndraws(res$parameters), 7)

  direct_func <- function(n_datasets, base_indices = 1:length(res)) {
      res[base_indices[rep(1:length(base_indices), length.out = n_datasets)]]
    }

  res_direct1 <- generate_datasets(SBC_generator_custom(direct_func),  n_datasets = 7)

  expect_equal(res, res_direct1, check.attributes = FALSE)

  res_direct2 <- generate_datasets(SBC_generator_custom(direct_func, base_indices = 1:3),
  n_datasets = 5)

  expect_identical(res[c(1,2,3,1,2)], res_direct2)

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
    SBC_generator_function(list_function, N = 10),
    n_datasets = 7)

  res_subs <- res[3:5]

  expect_identical(res_subs$parameters, posterior::subset_draws(res$parameters, draw = 3:5))
  expect_identical(res_subs$generated[[1]], res$generated[[3]])
  expect_identical(res_subs$generated[[2]], res$generated[[4]])
  expect_identical(res_subs$generated[[3]], res$generated[[5]])
  expect_equal(length(res_subs$generated), 3)

})
