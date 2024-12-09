test_that("combine_draws_matrix_for_bf", {
  dm0 <- posterior::draws_matrix("a" = c(1,2,3,4), "b" = c(5,6,7,8))
  dm1 <- posterior::draws_matrix("a" = c(10,20,30,40), "c" = c(50, 60, 70, 80))
  model_draws <- c(0, 1, 0, 1)

  res <- combine_draws_matrix_for_bf(dm0, dm1, model_draws, model_var = "model_test")
  target <- posterior::draws_matrix(
    model_test = model_draws,
    a = c(1,20,3,40),
    b = c(5, -Inf, 7, -Inf),
    c = c(-Inf, 60, -Inf, 80),
    .m0.a = dm0[,"a"],
    .m0.b = dm0[,"b"],
    .m1.a = dm1[,"a"],
    .m1.c = dm1[,"c"])
  expect_identical(res, target)

  res_NA_raw <- combine_draws_matrix_for_bf(dm0, dm1, model_draws, NA_raw_dm = TRUE, model_var = "model_test")
  target_NA_raw <- posterior::draws_matrix(
    model_test = model_draws,
    a = c(1,20,3,40),
    b = c(5, -Inf, 7, -Inf),
    c = c(-Inf, 60, -Inf, 80),
    .m0.a = c(1,NA,3,NA),
    .m0.b = c(5, NA, 7,NA),
    .m1.a = c(NA, 20, NA, 40),
    .m1.c = c(NA, 60, NA, 80))
  expect_identical(res_NA_raw, target_NA_raw)
})


test_that("combine_var_attributes_for_bf", {
  dm0 <- posterior::draws_matrix("a" = c(1,2,3,4), "b[0]" = c(5,6,7,8), "b[1]" = c(9,10,11,12))
  dm1 <- posterior::draws_matrix("a" = c(10,20,30,40), "c" = c(50, 60, 70, 80))
  var_attr0 <- var_attributes(a = binary_var_attribute())
  var_attr1 <- var_attributes(c = c(possibly_constant_var_attribute(), hidden_var_attribute()))

  expect_identical(
    combine_var_attributes_for_bf(dm0, dm1, var_attr0, var_attr1, model = "mmm"),
    var_attributes(
      .m0.a = c(hidden_var_attribute(), submodel_var_attribute(0), binary_var_attribute()),
      .m0.b = c(hidden_var_attribute(), submodel_var_attribute(0)),
      .m1.a = c(hidden_var_attribute(), submodel_var_attribute(1)),
      .m1.c = c(hidden_var_attribute(), submodel_var_attribute(1), possibly_constant_var_attribute(), hidden_var_attribute()),
      a = binary_var_attribute(),
      b = inf_valid_var_attribute(),
      c = c(inf_valid_var_attribute(), possibly_constant_var_attribute(), hidden_var_attribute()),
      mmm = c(binary_var_attribute(), possibly_constant_var_attribute())
    )
  )
})
