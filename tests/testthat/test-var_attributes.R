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
