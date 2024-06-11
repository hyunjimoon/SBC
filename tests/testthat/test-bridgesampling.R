test_that("combine_draws_matrix_for_bf", {
  dm0 <- posterior::draws_matrix("a" = c(1,2,3,4), "b" = c(5,6,7,8))
  dm1 <- posterior::draws_matrix("a" = c(10,20,30,40), "c" = c(50, 60, 70, 80))
  model_draws <- c(0, 1, 0, 1)

  res <- combine_draws_matrix_for_bf(dm0, dm1, model_draws, model_var = "model_test")
  target <- posterior::draws_matrix(
    model_test = model_draws,
    a = c(1,20,3,40),
    b = c(5, NA, 7, NA),
    c = c(NA, 60, NA, 80),
    .m0.a = dm0[,"a"],
    .m0.b = dm0[,"b"],
    .m1.a = dm1[,"a"],
    .m1.c = dm1[,"c"])
  expect_identical(res, target)

  res_NA_raw <- combine_draws_matrix_for_bf(dm0, dm1, model_draws, NA_raw_dm = TRUE, model_var = "model_test")
  target_NA_raw <- posterior::draws_matrix(
    model_test = model_draws,
    a = c(1,20,3,40),
    b = c(5, NA, 7, NA),
    c = c(NA, 60, NA, 80),
    .m0.a = c(1,NA,3,NA),
    .m0.b = c(5, NA, 7,NA),
    .m1.a = c(NA, 20, NA, 40),
    .m1.c = c(NA, 60, NA, 80))
  expect_identical(res_NA_raw, target_NA_raw)
})
