test_that("attribute_present", {
  expect_identical(
    attribute_present("binary",
                      c("a", "b"),
                      var_attributes = list(
                        b = c("allow_na", "binary"),
                        d = c("binary")
                      )),
    c("a" = FALSE, "b" = TRUE)
  )
  expect_identical(
    attribute_present("allow_na",
                      c("a", "b", "c", "d"),
                      var_attributes = list(
                        a = c("binary", "test"),
                        b = c("allow_na", "binary"),
                        d = c("binary"),
                        e = character(0)
                      )),
    c("a" = FALSE, "b" = TRUE, "c" = FALSE, "d" = FALSE)
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

  expect_false(attribute_present_stats("binary", "allow_na"))
  expect_false(attribute_present_stats("binary", "binary2"))
  expect_false(attribute_present_stats("binary", "allow_binary, other"))
  expect_false(attribute_present_stats("binary", "allow_binary, other"))
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

})
