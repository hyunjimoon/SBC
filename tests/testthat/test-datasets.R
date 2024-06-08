test_that("Generating datasets via functions", {
  list_function <- function(N) {
    mu <- rnorm(1)
    sigma <- abs(rnorm(1))
    beta <- matrix(1:6, nrow = 3, ncol = 2)
    gamma <- array(1:12, dim = c(2,3,2))
    delta <- rnorm(3)
    names(delta) <- LETTERS[1:3]
    y1 <- rnorm(N, mu, sigma)
    y2 <- rnorm(2 * N, mu + 5, sigma)
    list(variables = list(mu = mu, sigma = sigma, beta = beta, gamma = gamma, delta = delta),
         generated = list(y1 = y1, y2 = y2))
  }

  res <- generate_datasets(
    SBC_generator_function(list_function, N = 10),
    n_sims = 7)

  expect_true(length(res) == 7)

  beta_vars <- paste0("beta[", rep(1:3, times = 2), ",", rep(1:2, each = 3), "]")
  gamma_vars <- paste0("gamma[", rep(1:2, times = 6), ",", rep(rep(1:3, each = 2), times = 2),
                       ",", rep(1:2, each = 6), "]")
  delta_vars <- paste0("delta[",LETTERS[1:3],"]")
  expect_identical(posterior::variables(res$variables), c("mu", "sigma", beta_vars, gamma_vars, delta_vars))
  expect_identical(names(res$generated[[1]]), c("y1", "y2"))

  expect_equal(posterior::ndraws(res$variables), 7)

  direct_func <- function(n_sims, base_indices = 1:length(res)) {
      res[base_indices[rep(1:length(base_indices), length.out = n_sims)]]
    }

  res_direct1 <- generate_datasets(SBC_generator_custom(direct_func),  n_sims = 7)

  expect_equal(res, res_direct1, check.attributes = FALSE)

  res_direct2 <- generate_datasets(SBC_generator_custom(direct_func, base_indices = 1:3),
  n_sims = 5)

  expect_identical(res[c(1,2,3,1,2)], res_direct2)

})

test_that("Generating datasets via functions - exceptions", {
  missing_gen_function <- function() {
    list(variables = list(mu = 1),
         not_generated = 1)
  }

  expect_error(generate_datasets(
    SBC_generator_function(missing_gen_function),
    n_sims = 1), class = "SBC_datasets_error")

  missing_par_function <- function() {
    list(not_variables = list(mu = 1),
         generated = 1)
  }

  expect_error(generate_datasets(
    SBC_generator_function(missing_par_function),
    n_sims = 1), class = "SBC_datasets_error")


  missing_names_function <- function() {
    list(variables = list(mu = 1, 5),
         generated = 1)
  }

  expect_error(generate_datasets(
    SBC_generator_function(missing_names_function),
    n_sims = 1), class = "SBC_datasets_error")
})

test_that("subsetting datasets", {
  my_var_attributes <- list(
    mu = hidden_var_attribute()
  )
  list_function <- function(N) {
    mu <- rnorm(1)
    sigma <- abs(rnorm(1))
    y1 <- rnorm(N, mu, sigma)
    y2 <- rnorm(2 * N, mu + 5, sigma)
    list(variables = list(mu = mu, sigma = sigma),
         generated = list(y1 = y1, y2 = y2),
         var_attributes = my_var_attributes)
  }

  res <- generate_datasets(
    SBC_generator_function(list_function, N = 10),
    n_sims = 7)

  res_subs <- res[3:5]

  expect_identical(res_subs$variables, posterior::subset_draws(res$variables, draw = 3:5))
  expect_identical(res_subs$generated[[1]], res$generated[[3]])
  expect_identical(res_subs$generated[[2]], res$generated[[4]])
  expect_identical(res_subs$generated[[3]], res$generated[[5]])
  expect_identical(res_subs$var_attributes, res$var_attributes)
  expect_equal(length(res_subs$generated), 3)

})

test_that("binding datasets", {
  my_var_attributes <- list(
    mu = hidden_var_attribute()
  )
  list_function <- function(N) {
    mu <- rnorm(1)
    sigma <- abs(rnorm(1))
    y1 <- rnorm(N, mu, sigma)
    y2 <- rnorm(2 * N, mu + 5, sigma)
    list(variables = list(mu = mu, sigma = sigma),
         generated = list(y1 = y1, y2 = y2),
         var_attributes = my_var_attributes)
  }

  res <- generate_datasets(
    SBC_generator_function(list_function, N = 10),
    n_sims = 7)

  res2 <- generate_datasets(
    SBC_generator_function(list_function, N = 10),
    n_sims = 5)


  res_combined <- bind_datasets(res, res2)

  expect_identical(res$var_attributes, res_combined$var_attributes)
  expect_equal(length(res_combined$generated), 12)

  res_no_attributes <- res2
  res_no_attributes$var_attributes <- NULL
  expect_error(bind_datasets(res, res_no_attributes), "var_attributes")

  res_other_attributes <- res2
  res_other_attributes$var_attributes <- list(sigma = possibly_constant_var_attribute())
  expect_error(bind_datasets(res, res_other_attributes), "var_attributes")

})
