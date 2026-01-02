test_that("attribute_present", {
  expect_identical(
    attribute_present("binary",
                      c("a", "b"),
                      var_attr = var_attributes(
                        b = c("allow_na", "binary"),
                        d = c("binary")
                      )),
    c("a" = FALSE, "b" = TRUE)
  )
  expect_identical(
    attribute_present("allow_na",
                      c("a", "b", "b[1]", "b[2,3,4]", "c", "d"),
                      var_attr = var_attributes(
                        a = c("binary", "test"),
                        b = c("allow_na", "binary"),
                        d = c("binary"),
                        e = character(0)
                      )),
    c("a" = FALSE, "b" = TRUE, "b[1]" = TRUE, "b[2,3,4]" = TRUE,"c" = FALSE, "d" = FALSE)
  )
})

test_that("attribute_present_stats", {
  expect_true(attribute_present_stats("binary", "binary"))
  expect_true(attribute_present_stats("binary", "binary, allow_na"))
  expect_true(attribute_present_stats("binary", "binary,allow_na"))
  expect_true(attribute_present_stats("binary", "allow_na,binary"))
  expect_true(attribute_present_stats("binary", "allow_na, binary"))
  expect_true(attribute_present_stats("binary", "allow_na, binary, other"))
  expect_true(attribute_present_stats("binary", "allow_na,binary,other"))
  expect_true(attribute_present_stats("binary", "allow_na,binary, other"))
  expect_true(attribute_present_stats("binary", "allow_na, binary, other"))
  expect_true(attribute_present_stats("binary", "binary, binary, binary"))
  expect_true(attribute_present_stats("submodel(0)", "binary, binary, submodel(0)"))
  expect_true(attribute_present_stats("submodel(1)", "submodel(1), binary, binary"))

  expect_false(attribute_present_stats("binary", "allow_na"))
  expect_false(attribute_present_stats("binary", "binary2"))
  expect_false(attribute_present_stats("binary", "allow_binary, other"))
  expect_false(attribute_present_stats("binary", "allow_binary, other"))
  expect_false(attribute_present_stats("submodel(0)", "submodel(1), binary, binary"))
})


test_that("extract_attribute_arguments_stats", {
  expect_identical(extract_attribute_arguments_stats("submodel", c("binary", "submodel", "submodel(1)", "submodel(2,3,4)")),
                   c(NA_character_, NA_character_, "1", "2,3,4"))
  expect_identical(extract_attribute_arguments_stats("submodel", c("binary,allow_na", "notsubmodel(1)", "a,submodel(1)", "binary, submodel(2,3,4)")),
                   c(NA_character_, NA_character_, "1", "2,3,4"))
  expect_identical(extract_attribute_arguments_stats("submodel", c("binary", "submodel,questor", "aa,submodel(1), binary(3)", "qqq(a,b),submodel(2,3,4), prot(3,2)")),
                   c(NA_character_, NA_character_, "1", "2,3,4"))
  expect_identical(extract_attribute_arguments_stats("submodel", c("binary", "submodel,questor", "aa,submodel(1), binary(3)", "qqq(a,b),submodel(3), prot(3,2),submodel(2,3,4)")),
                   c(NA_character_, NA_character_, "1", "2,3,4"))
})



test_that("var_attributes_to_attributes_column", {
  expect_identical(
    var_attributes_to_attributes_column(list(a = "discrete"), c("a", "b")),
    c("discrete", "")
  )

  expect_identical(
    var_attributes_to_attributes_column(list(a = "discrete", b = c("test1", "test2", "test3")), c("c", "d", "b", "a")),
    c("", "", "test1, test2, test3", "discrete")
  )

  expect_identical(
    var_attributes_to_attributes_column(list(a = "discrete", b = c("test1", "test2", "test3")), c("a", "ab", "bb", "a")),
    c("discrete", "", "", "discrete")
  )

  expect_identical(
    var_attributes_to_attributes_column(NULL, c("a", "ab", "bb", "a")),
    c("", "", "", "")
  )

  expect_identical(
    var_attributes_to_attributes_column(list(a = "discrete", b = c("test1", "test2", "test3")), c("a[1]", "a[2]", "b[1,1]", "b[4,5]")),
    c("discrete", "discrete", "test1, test2, test3", "test1, test2, test3")
  )

})

test_that("combine_var_attributes", {
  expect_identical(combine_var_attributes(
    var_attributes(a = c(binary_var_attribute(), hidden_var_attribute()),
                   b = c(hidden_var_attribute())),
    var_attributes(a = c(hidden_var_attribute()),
                   q = possibly_constant_var_attribute(),
                   z = binary_var_attribute()),
    var_attributes(a  = hidden_var_attribute(),
                   b = binary_var_attribute(),
                   q = na_valid_var_attribute())
  ),
  var_attributes(a = c(binary_var_attribute(), hidden_var_attribute(),hidden_var_attribute(),hidden_var_attribute()),
                 b = c(hidden_var_attribute(), binary_var_attribute()),
                 q = c(possibly_constant_var_attribute(), na_valid_var_attribute()),
                 z = binary_var_attribute()))
})

test_that("remove_attribute_from_stats", {
  expect_identical(remove_attribute_from_stats("a", "a, b,c"), "b,c")
  expect_identical(remove_attribute_from_stats("b", "ab,b, cb,b,bbc"), "ab, cb,b,bbc")
  expect_identical(remove_attribute_from_stats("ab", "a,b,cab,ab"), "a,b,cab")
  expect_identical(remove_attribute_from_stats("ab(44)", "a,b,cab,ab(44)"), "a,b,cab")
  expect_identical(remove_attribute_from_stats("ab(44)", "a,b,cab,ab(21)"), "a,b,cab,ab(21)")
})


test_that("compute_default_diagnostics_and_attributes", {
  example_fit <- posterior::draws_matrix(only_inf = rep(Inf, 300),
                                         one_inf = c(Inf, rnorm(299)),
                                         only_NA = rep(NA, 300),
                                         one_NA = c(NA, rnorm(299)))

  variables <- posterior::draws_matrix(only_inf = Inf,
                                       one_inf = 0,
                                       only_NA = NA,
                                       one_NA = 0)


  stats_no_attr <-
    SBC_statistics_from_single_fit(example_fit, variables = variables,
                                 generated = list(),
                                 thin_ranks = 1,
                                 ensure_num_ranks_divisor = 1,
                                 dquants = NULL,
                                 backend = NULL
                                 )

  expect_identical(is.na(stats_no_attr$q5), c(F, F, T, T))
  expect_identical(is.na(stats_no_attr$rhat), c(T, F, T, T))
  expect_identical(is.na(stats_no_attr$ess_bulk), c(T, F, T, T))
  expect_identical(is.na(stats_no_attr$ess_tail), c(T, T, T, T))

  expect_identical(stats_no_attr$has_na, c(F, F, T, T))
  expect_identical(stats_no_attr$all_inf, c(T, F, F, F))

  stats_attr <-
    SBC_statistics_from_single_fit(example_fit, variables = variables,
                                 generated = list(),
                                 thin_ranks = 1,
                                 ensure_num_ranks_divisor = 1,
                                 dquants = NULL,
                                 backend = NULL,
                                 var_attributes = var_attributes(
                                   only_inf = inf_valid_var_attribute(),
                                   one_inf = inf_valid_var_attribute(),
                                   only_NA = na_valid_var_attribute(),
                                   one_NA = na_valid_var_attribute())
  )


  expect_identical(is.na(stats_attr$q5), c(F, F, T, F))
  expect_identical(is.na(stats_attr$rhat), c(T, F, T, F))
  expect_identical(is.na(stats_attr$ess_bulk), c(T, F, T, F))
  expect_identical(is.na(stats_attr$ess_tail), c(T, F, T, F))

  expect_identical(stats_attr$has_na, c(F, F, T, T))
  expect_identical(stats_attr$all_inf, c(T, F, F, F))


  stats_no_attr$sim_id <- 1
  diag_no_attr <- compute_default_diagnostics(stats_no_attr)
  expect_equal(diag_no_attr$n_has_na, 2)
  expect_equal(diag_no_attr$n_na_rhat, 3)
  expect_equal(diag_no_attr$n_na_ess_bulk, 3)
  expect_equal(diag_no_attr$n_na_ess_tail, 4)
  expect_true(is.na(diag_no_attr$min_ess_bulk))
  expect_true(is.na(diag_no_attr$min_ess_tail))
  expect_true(is.na(diag_no_attr$min_ess_to_rank))

  stats_attr$sim_id <- 1
  diag_attr <- compute_default_diagnostics(stats_attr)
  expect_equal(diag_attr$n_has_na, 0)
  expect_equal(diag_attr$n_na_rhat, 0)
  expect_equal(diag_attr$n_na_ess_bulk, 0)
  expect_equal(diag_attr$n_na_ess_tail, 0)
  expect_false(is.na(diag_attr$min_ess_bulk))
  expect_false(is.na(diag_attr$min_ess_tail))
  expect_false(is.na(diag_attr$min_ess_to_rank))


})
