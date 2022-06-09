test_that("empirical coverage and merging", {
  stats_rnd <- tidyr::crossing(variable = LETTERS[1:5], max_rank = 100, sim_id = 1:100)
  stats_rnd <- dplyr::mutate(stats_rnd, rank = sample(0:100, size = dplyr::n(), replace = TRUE, prob = c(1:50,seq(200, 40, length.out = 51))))
  #widths_to_test <- c(0.3,0.6,0.9)
  widths_to_test <- c(0.3)
  per_var <- empirical_coverage(stats_rnd, width = widths_to_test)
  stats_merged <- dplyr::mutate(stats_rnd, variable = "Merged")
  merged <- empirical_coverage(stats_merged, width = widths_to_test)

  per_var_merged <- dplyr::summarise(
    dplyr::group_by(per_var, width),
    estimate_mean = mean(estimate),
    min_ci_low = min(ci_low),
    max_ci_high = max(ci_high)
  )

  expect_equal(merged$estimate, per_var_merged$estimate_mean)
  for(i in length(widths_to_test)) {
    expect_gt(merged$ci_low[i], per_var_merged$min_ci_low[i])
    expect_lt(merged$ci_high[i], per_var_merged$max_ci_high[i])
  }
})
