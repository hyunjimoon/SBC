test_that("Response sequences for univariate and multivariate models", {
  # TODO: the order within each list element shouldn't matter,
  # but now tests are sensitive to change in that order.
  # Keep that in mind when investigating test failures

  expect_equal(brms_response_sequence(
    brms::bf(y ~ x)),
    list("y"))

  expect_equal(brms_response_sequence(
    brms::bf(y ~ x) + brms::bf(y2 ~ x) + brms::set_rescor(FALSE)),
    list(c("y", "y2")))

  expect_equal(brms_response_sequence(
    brms::bf(z ~ y) + brms::bf(y ~ x) + brms::set_rescor(FALSE)),
    list("y", "z"))

  expect_equal(brms_response_sequence(
    brms::bf(q ~ y, shape ~ z, family = "negbinomial") + brms::bf(z ~ y + r) + brms::bf(y ~ x) + brms::bf(r ~ x) + brms::set_rescor(FALSE)),
    list(c("y", "r"), "z", "q"))

  expect_equal(brms_response_sequence(
    brms::bf(q ~ y, shape ~ x, family = "negbinomial") + brms::bf(z ~ y + r) + brms::bf(y ~ x) + brms::bf(r ~ x) + brms::set_rescor(FALSE)),
    list(c("y", "r"), c("q", "z")))
})
