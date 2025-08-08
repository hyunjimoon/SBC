test_that("SBC_backend_cached", {
  base_backend <- SBC_backend_mock(output = "OUT34885", warning = "WARN24785", bridgesampler = list(5))

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
  expect_true(file.exists(cache_file_test))
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


    expect_message(
      regexp = "Storing bridge",
      SBC_fit_to_bridge_sampler(cached_backend, fit = list(), generated = list())
      )


  cache_file_test_bs <- cached_fit_filename(cache_dir, base_backend, list(), suffix = "_bridgesampling")
  expect_true(file.exists(cache_file_test_bs))
  cached_bs <- readRDS(cache_file_test_bs)
  expect_identical(cached_bs$bridgesampler, base_backend$bridgesampler)

    expect_message(
      regexp = "bridge.*read from cache",
      SBC_fit_to_bridge_sampler(cached_backend, fit = list(), generated = list())
      )
})
