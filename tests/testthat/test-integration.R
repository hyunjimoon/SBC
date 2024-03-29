test_that("Integration test with mock backend", {
  set.seed(546622)
  N_sims <- 300
  a_vals <- rnorm(N_sims)
  b_vals <- rgamma(N_sims, shape = 3)
  fit_result <- posterior::draws_matrix(a = a_vals, b = b_vals)
  backend <- SBC_backend_mock(result = fit_result,
                              output = "OUT",
                              message = "MSG",
                              warning = "WARN"
                              )

  # Reuse `fit_result` + shift them slightly as the true values, ensuring uniform ranks
  a_vals_true <- a_vals + min(abs(diff(sort(a_vals)))) / 2
  b_vals_true <- b_vals + min(abs(diff(sort(b_vals)))) / 2
  true_result <- posterior::draws_matrix(a = a_vals_true, b = b_vals_true)
  datasets <- SBC_datasets(true_result,generated = rep(list(NULL), N_sims))

  res <- compute_SBC(datasets, backend, thin_ranks = 1,
                     ensure_num_ranks_divisor = 1)

  expected_ranks <- rep(1:N_sims, each = 2)
  expect_equivalent(sort(res$stats$rank), expected_ranks)
  expect_identical(res$outputs, rep(list("OUT"), N_sims))
  expect_identical(res$messages, rep(list("MSG\n"), N_sims))
  expect_identical(res$warnings, rep(list("WARN"), N_sims))
  expect_identical(res$errors, rep(list(NULL), N_sims))
  expect_identical(res$fits, rep(list(fit_result), N_sims))

  backend2 <- backend
  backend2$error <- SBC:::SBC_error("SBC_test_error", "ERR")
  res2_with_outputs <- SBC:::capture_all_outputs(compute_SBC(datasets, backend2, thin_ranks = 1))
  res2 <- res2_with_outputs$result
  expect_identical(res2$errors, rep(list(backend2$error), N_sims))
  expect_identical(res2$fits, rep(list(NULL), N_sims))

  expect_equal(res2_with_outputs$warnings, "All simulations produced error when fitting")

})


test_that("Result caching", {
  set.seed(1521336)
  N_sims <- 100
  a_vals <- rnorm(N_sims)
  fit_result <- posterior::draws_matrix(a = a_vals)

  backend <- SBC_backend_mock(result = fit_result)
  datasets <- SBC_datasets(fit_result, generated = rep(list(NULL), N_sims))


  cache_file <- tempfile(fileext = ".rds")

  res_first <- SBC:::capture_all_outputs(
    compute_SBC(datasets, backend, thin_ranks = 1, cache_mode = "results", cache_location = cache_file))

  expect_false(any(grepl("cache",  c(res_first$output, res_first$messages, res_first$warnings))))

  file_info_cols <- c("size", "mtime", "ctime")

  # Succesful load from cache
  old_file_info <- file.info(cache_file)
  expect_message(
    compute_SBC(datasets, backend, thin_ranks = 1, cache_mode = "results", cache_location = cache_file),
    "loaded from cache"
    )
  new_file_info <- file.info(cache_file)
  expect_identical(old_file_info[, file_info_cols], new_file_info[, file_info_cols])

  # Change datasets
  old_file_info <- file.info(cache_file)
  datasets_changed <- datasets
  datasets_changed[[3]] <- "a"
  expect_message(
    compute_SBC(datasets_changed, backend, thin_ranks = 1, cache_mode = "results", cache_location = cache_file),
    "datasets.*differ.*recompute"
    )
  new_file_info <- file.info(cache_file)
  expect_lt(old_file_info$mtime, new_file_info$mtime)


  # Now should be succesful
  old_file_info <- file.info(cache_file)
  expect_message(
    compute_SBC(datasets_changed, backend, thin_ranks = 1, cache_mode = "results", cache_location = cache_file),
    "loaded from cache"
  )
  new_file_info <- file.info(cache_file)
  expect_identical(old_file_info[, file_info_cols], new_file_info[, file_info_cols])

  # Change backend
  old_file_info <- file.info(cache_file)
  backend_changed <- backend
  backend_changed$result[5, "a"] <- 0
  expect_message(
    compute_SBC(datasets_changed, backend_changed, thin_ranks = 1, cache_mode = "results", cache_location = cache_file),
    "backend.*differ.*recompute"
  )
  new_file_info <- file.info(cache_file)
  expect_lt(old_file_info$mtime, new_file_info$mtime)

  # Change thin_ranks, so should call recompute_statistics and not touch the cache
  old_file_info <- file.info(cache_file)
  expect_message(
    compute_SBC(datasets_changed, backend_changed, thin_ranks = 2, cache_mode = "results", cache_location = cache_file),
    "recompute_SBC_statistics"
  )
  new_file_info <- file.info(cache_file)
  expect_identical(old_file_info[, file_info_cols], new_file_info[, file_info_cols])

  # Repeating should give the same result
  expect_message(
    compute_SBC(datasets_changed, backend_changed, thin_ranks = 2, cache_mode = "results", cache_location = cache_file),
    "recompute_SBC_statistics"
  )
  new_file_info <- file.info(cache_file)
  expect_identical(old_file_info[, file_info_cols], new_file_info[, file_info_cols])

  # Start fresh again
  # Store cache but with keep_fits = FALSE
  cache_file2 <- tempfile(fileext = ".rds")

  res_first_nokeep <- SBC:::capture_all_outputs(
    compute_SBC(datasets, backend, thin_ranks = 1, cache_mode = "results", keep_fits = FALSE, cache_location = cache_file2))

  expect_false(any(grepl("cache",  c(res_first_nokeep$output, res_first_nokeep$messages, res_first_nokeep$warnings))))

  # Change thin_ranks, but since keep_fits = FALSE should force a refit and update cache
  old_file_info <- file.info(cache_file2)
  expect_message(
    compute_SBC(datasets, backend, thin_ranks = 2, cache_mode = "results", cache_location = cache_file2),
    "exists.*keep_fits == FALSE.*recompute"
  )
  new_file_info <- file.info(cache_file2)
  expect_lt(old_file_info$mtime, new_file_info$mtime)

  # Repeating should give the same result
  old_file_info <- file.info(cache_file2)
  expect_message(
    compute_SBC(datasets, backend, thin_ranks = 2, cache_mode = "results", cache_location = cache_file2),
    "loaded from cache"
  )
  new_file_info <- file.info(cache_file2)
  expect_identical(old_file_info[, file_info_cols], new_file_info[, file_info_cols])

})
