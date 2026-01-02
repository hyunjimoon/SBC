test_that("combine_draws_matrix_for_bf", {
  dm0 <- posterior::draws_matrix("a" = c(1,2,3,4), "b" = c(5,6,7,8))
  dm1 <- posterior::draws_matrix("a" = c(10,20,30,40), "c" = c(50, 60, 70, 80))
  model_draws <- c(0, 1, 0, 1)

  res <- combine_draws_matrix_for_bf(list(dm0, dm1), model_draws, model_var = "model_test")
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

  res_NA_raw <- combine_draws_matrix_for_bf(list(dm0, dm1), model_draws, NA_raw_dm = TRUE, model_var = "model_test")
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


  dm2 <- posterior::draws_matrix("a" = c(100,200,300,400), "b" = c(500, 600, 700, 800), "d" = c(900, 1000, 1100, 1200))
  model_draws3 <- c(1, 2, 0, 2)

  res3 <- combine_draws_matrix_for_bf(list(dm0, dm1, dm2), model_draws3, model_var = "model_test")
  target3 <- posterior::draws_matrix(
    model_test = model_draws3,
    is_model_test0 = c(0, 0, 1, 0),
    is_model_test1 = c(1, 0, 0, 0),
    is_model_test2 = c(0, 1, 0, 1),
    a = c(10,200,3,400),
    b = c(-Inf, 600, 7, 800),
    c = c(50, -Inf, -Inf, -Inf),
    d = c(-Inf, 1000, -Inf, 1200),
    .m0.a = dm0[,"a"],
    .m0.b = dm0[,"b"],
    .m1.a = dm1[,"a"],
    .m1.c = dm1[,"c"],
    .m2.a = dm2[,"a"],
    .m2.b = dm2[,"b"],
    .m2.d = dm2[,"d"]
  )
  expect_identical(res3, target3)

  res3_NA_raw <- combine_draws_matrix_for_bf(list(dm0, dm1, dm2), model_draws3, NA_raw_dm = TRUE, model_var = "model_test")
  target3_NA_raw <- posterior::draws_matrix(
    model_test = model_draws3,
    is_model_test0 = c(0, 0, 1, 0),
    is_model_test1 = c(1, 0, 0, 0),
    is_model_test2 = c(0, 1, 0, 1),
    a = c(10,200,3,400),
    b = c(-Inf, 600, 7, 800),
    c = c(50, -Inf, -Inf, -Inf),
    d = c(-Inf, 1000, -Inf, 1200),
    .m0.a = if_else(model_draws3 == 0, dm0[,"a"], NA),
    .m0.b = if_else(model_draws3 == 0, dm0[,"b"], NA),
    .m1.a = if_else(model_draws3 == 1, dm1[,"a"], NA),
    .m1.c = if_else(model_draws3 == 1, dm1[,"c"], NA),
    .m2.a = if_else(model_draws3 == 2, dm2[,"a"], NA),
    .m2.b = if_else(model_draws3 == 2, dm2[,"b"], NA),
    .m2.d = if_else(model_draws3 == 2, dm2[,"d"], NA))
  expect_identical(res3_NA_raw, target3_NA_raw)
})


test_that("combine_var_attributes_for_bf", {
  dm0 <- posterior::draws_matrix("a" = c(1,2,3,4), "b[0]" = c(5,6,7,8), "b[1]" = c(9,10,11,12))
  dm1 <- posterior::draws_matrix("a" = c(10,20,30,40), "c" = c(50, 60, 70, 80))

  var_attr0 <- var_attributes(a = binary_var_attribute())
  var_attr1 <- var_attributes(c = c(possibly_constant_var_attribute(), hidden_var_attribute()))

  expect_identical(
    combine_var_attributes_for_bf(list(dm0, dm1), list(var_attr0, var_attr1), model = "mmm"),
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

  dm2 <- posterior::draws_matrix("a" = c(100,200,300,400), "b" = c(500, 600, 700, 800), "d" = c(900, 1000, 1100, 1200))
  var_attr2 <- var_attributes(a = inf_valid_var_attribute(), d = c(possibly_constant_var_attribute(), hidden_var_attribute()))

  expect_identical(
    combine_var_attributes_for_bf(list(dm0, dm1, dm2), list(var_attr0, var_attr1, var_attr2), model = "mmm"),
    var_attributes(
      .m0.a = c(hidden_var_attribute(), submodel_var_attribute(0), binary_var_attribute()),
      .m0.b = c(hidden_var_attribute(), submodel_var_attribute(0)),
      .m1.a = c(hidden_var_attribute(), submodel_var_attribute(1)),
      .m1.c = c(hidden_var_attribute(), submodel_var_attribute(1), possibly_constant_var_attribute(), hidden_var_attribute()),
      .m2.a = c(hidden_var_attribute(), submodel_var_attribute(2), inf_valid_var_attribute()),
      .m2.b = c(hidden_var_attribute(), submodel_var_attribute(2)),
      .m2.d = c(hidden_var_attribute(), submodel_var_attribute(2), possibly_constant_var_attribute(), hidden_var_attribute()),
      a = c(binary_var_attribute(), inf_valid_var_attribute()),
      b = inf_valid_var_attribute(),
      c = c(inf_valid_var_attribute(), possibly_constant_var_attribute(), hidden_var_attribute()),
      d = c(inf_valid_var_attribute(), possibly_constant_var_attribute(), hidden_var_attribute()),
      is_mmm0 = c(binary_var_attribute(), possibly_constant_var_attribute()),
      is_mmm1 = c(binary_var_attribute(), possibly_constant_var_attribute()),
      is_mmm2 = c(binary_var_attribute(), possibly_constant_var_attribute()),
      mmm = c(possibly_constant_var_attribute())
    )
  )

})


test_that("bridgesampling_diagnostics_special_treatment", {
  diags1 <- structure(data.frame(sim_id = 1L, prob_H1 = 0.5, test_H0 = 0.1, test_H1 = 0.5), class = c("SBC_bridgesampling_diagnostics", "data.frame"))
  attr(diags1, "submodel_classes") <- list("H0" = c("test_class", "data.frame"), "H1" = c("test_class2", "data.frame"))
  diags2 <- structure(data.frame(sim_id = 2L, prob_H1 = 0.3, test_H0 = 0.3, test_H1 = 0.3), class = c("SBC_bridgesampling_diagnostics", "data.frame"))
  attr(diags2, "submodel_classes") <- list("H0" = c("test_class", "data.frame"), "H1" = c("test_class2", "data.frame"))
  diags3 <- structure(data.frame(sim_id = 3L, prob_H1 = 0.1, test_H0 = 0.5, test_H1 = 0.1), class = c("SBC_bridgesampling_diagnostics", "data.frame"))
  attr(diags3, "submodel_classes") <- list("H0" = c("test_class", "data.frame"), "H1" = c("test_class2", "data.frame"))

  diags_bound <- rbind(diags1, diags2, diags3)

  diags_expected <- structure(data.frame(sim_id = 1:3, prob_H1 = c(0.5, 0.3,0.1), test_H0 = c(0.1, 0.3, 0.5), test_H1 = c(0.5, 0.3, 0.1)), class = c("SBC_bridgesampling_diagnostics", "data.frame"))
  attr(diags_expected, "submodel_classes") <- list("H0" = c("test_class", "data.frame"), "H1" = c("test_class2", "data.frame"))

  expect_identical(diags_bound, diags_expected)

  diags_mismatched <- structure(data.frame(sim_id = 4L, prob_H1 = 1, test_H0 = 1, test_H1 = 1), class = c("SBC_bridgesampling_diagnostics", "data.frame"))
  attr(diags_mismatched, "submodel_classes") <- list("H0" = c("other_class", "data.frame"), "H1" = c("test_class2", "data.frame"))

  diags_bound_mismatch <- expect_warning(rbind(diags1, diags2, diags3, diags_mismatched), "Non-unique submodel classes")

  diags_expected_mismatch <- structure(data.frame(sim_id = 1:4, prob_H1 = c(0.5, 0.3,0.1, 1), test_H0 = c(0.1, 0.3, 0.5, 1), test_H1 = c(0.5, 0.3, 0.1, 1)), class = c("SBC_bridgesampling_diagnostics", "data.frame"))
  attr(diags_expected_mismatch, "submodel_classes") <- list("H0" = c("test_class", "data.frame"), "H1" = c("test_class2", "data.frame"))

  expect_identical(diags_bound_mismatch, diags_expected_mismatch)

  library(dplyr)
  diags_selected <- select(diags_bound, -sim_id)
  diags_selected_expected <- select(as.data.frame(diags_bound), -sim_id)
  class(diags_selected_expected) <- c("SBC_bridgesampling_diagnostics", "data.frame")
  attr(diags_selected_expected, "submodel_classes") <- list("H0" = c("test_class", "data.frame"), "H1" = c("test_class2", "data.frame"))

  expect_identical(diags_selected, diags_selected_expected)

})

test_that("discrete_to_cdf", {
  expect_equal(
    binary_to_cdf("test2", 0.4, 1),
    data.frame(variable = "test2", cdf_low = 0.6, cdf_high = 1)
  )
  expect_equal(
    binary_to_cdf("test", 0.65, 0),
    data.frame(variable = "test", cdf_low = 0, cdf_high = 0.35)
  )
  expect_equal(
    discrete_to_cdf("test", c(0.1, 0.2, 0.3, 0.2, 0.2), 4),
    data.frame(variable = "test", cdf_low = 0.8 , cdf_high = 1)
  )
  expect_equal(
    discrete_to_cdf("test", c(0.1, 0.2, 0.3, 0.2, 0.2), 2),
    data.frame(variable = "test", cdf_low = 0.1 + 0.2 , cdf_high = 0.1+0.2+0.3)
  )
  expect_equal(
    discrete_to_cdf("test", c(0.1, 0.2, 0.3, 0.2, 0.2), 0),
    data.frame(variable = "test", cdf_low = 0.0 , cdf_high = 0.1)
  )

})
