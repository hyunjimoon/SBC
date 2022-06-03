test_that("bind_generated_quantities", {
  expect_identical(
    bind_generated_quantities(
      generated_quantities(a = 1, .globals = "c"),
      generated_quantities(b = a + sqrt(y), .globals = list("x" = function(y) y ^ 2))
    ),
    generated_quantities(a = 1, b = a + sqrt(y), .globals = list("c", "x" = function(y) y ^ 2))
  )
})
