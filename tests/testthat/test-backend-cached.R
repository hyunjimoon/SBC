test_that("SBC_backend_cached", {
  base_backend <- SBC_backend_mock(output = "OUT34885", warning = "WARN24785")
  base_backend$result

  cache_dir <- tempfile()
  dir.create(cache_dir)
  cached_backend <- SBC_backend_cached(cache_dir, base_backend)

  expect_warning(
    regexp = "WARN24785",
    expect_message(
      regexp = "Storing fit",
      expect_output(
        regexp = c("OUT34885"),
        SBC_fit(cached_backend, list(), 1)
  )))


  cache_file_test <- cached_fit_filename(cache_dir, base_backend, list())
  expect_true(file.exists(
    cached_fit_filename(cache_dir, base_backend, list())
  ))
  cc <- readRDS(cache_file_test)
  expect_identical(cc$fit, base_backend$result)

  expect_warning(
    regexp = "WARN24785",
    expect_message(
      regexp = "read from cache",
      expect_output(
        regexp = c("OUT34885"),
        SBC_fit(cached_backend, list(), 1)
      )))

})
