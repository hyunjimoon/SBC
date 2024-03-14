test_that("combine_XXX functions work", {
  var_names_to_test <- c("A", "B[1]", "B[2]", "C[1,1]", "C[2,2]", "C[3,3]")
  expect_identical(combine_all_variables(var_names_to_test),
                   list(all = var_names_to_test))

  expect_identical(combine_array_elements(var_names_to_test),
                   list(A = "A", "B[]" = c("B[1]", "B[2]"), "C[]" = c("C[1,1]", "C[2,2]", "C[3,3]")))

})
