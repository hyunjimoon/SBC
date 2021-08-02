test_that("combine_args", {
  expect_identical(combine_args(list(a = 1, b = 2, c = 3),
                                list(b = "no", d = 4)),
                   list(a = 1, b = "no", c = 3, d = 4))

  expect_identical(combine_args(list(1, 2),
                                list(c = 3, d = 4)),
                   list(1, 2, c = 3, d = 4))

  expect_identical(combine_args(list(1, b = 2),
                                list(c = 3, d = 4)),
                   list(1, b = 2, c = 3, d = 4))

  expect_identical(combine_args(list(1, b = 2, c = 3),
                                list(13, b = "ugh", e = 5)),
                   list(1, b = "ugh", c = 3, 13, e = 5))

  expect_identical(combine_args(list(a = 1, b = 2, c = 3),
                                list(c = "ugh", a = "no")),
                   list(a = "no", b = 2, c = "ugh"))

})
