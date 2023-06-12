test_that("Response sequences for univariate and multivariate models", {
  # TODO: the order within each list element shouldn't matter,
  # but now tests are sensitive to change in that order.
  # Keep that in mind when investigating test failures
  expect_equal(brms_response_sequence(
    bf(z ~ y) + bf(y ~ x) + set_rescor(FALSE)),
    list("y", "z"))

  expect_equal(brms_response_sequence(
    bf(q ~ y, shape ~ z, family = "negbinomial") + bf(z ~ y + r) + bf(y ~ x) + bf(r ~ x) + set_rescor(FALSE)),
    list(c("y", "r"), "z", "q"))

  expect_equal(brms_response_sequence(
    bf(q ~ y, shape ~ x, family = "negbinomial") + bf(z ~ y + r) + bf(y ~ x) + bf(r ~ x) + set_rescor(FALSE)),
    list(c("y", "r"), c("q", "z")))
})
