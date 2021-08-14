test_that("capture_all_outputs", {
    expect_identical(
        capture_all_outputs({
            cat("Test")
            warning("W")
            message("M")
            14
            }),
        list(result = 14,
             messages = data.frame(type = c("warning", "message"),
                                   message = c("W", "M\n")),
             output = "Test"))
})

test_that("calculate_ranks_draws_matrix works", {

    dm <- matrix(NA_real_, nrow = 10, ncol = 4)
    colnames(dm) <- c("a","b","c", "d")

    dm[, "a"] <- sample(1:10)
    dm[, "b"] <- sample(1:10)
    dm[, "c"] <- sample(1:10)
    dm[, "d"] <- sample(c(1:5, 1:5))
    dm <- posterior::as_draws_matrix(dm)

    params <- matrix(c(3.5, -5, 15, 3), nrow = 1)
    colnames(params) <- c("a","b","c", "d")

    N_steps <- 1e4
    all_ranks <- matrix(NA_real_, nrow = N_steps, ncol = 4)
    for(i in 1:N_steps) {
        last_ranks <- calculate_ranks_draws_matrix(params, dm)
        all_ranks[i,] <- last_ranks

    }
    expect_true(all(all_ranks[,1] == 3))
    expect_true(all(all_ranks[,2] == 0))
    expect_true(all(all_ranks[,3] == 10))

    # The final rank is stochastic due to presence of ties
    expect_true(all(all_ranks[,4] >= 4))
    expect_true(all(all_ranks[,4] <= 6))

    rank4_stats <- table(all_ranks[,4])
    expect_true(all(rank4_stats > 0))
    expect_gt(chisq.test(rank4_stats)$p.val, 1e-15)

    expect_equal(length(last_ranks), 4)
    expect_equal(attr(last_ranks, "max_rank"), 10)

    expect_equal(names(last_ranks), posterior::variables(dm))
})

test_that("calculate_sds_draws_matrix", {
    dm <- matrix(NA_real_, nrow = 10, ncol = 3)
    colnames(dm) <- c("a","b","c")

    dm[, "a"] <- sample(1:10)
    dm[, "b"] <- sample(-4:5)
    dm[, "c"] <- sample(11:20)
    dm <- posterior::as_draws_matrix(dm)

    expected_res <- c(a = sd(1:10), b = sd(-4:5), c = sd(11:20))

    expect_identical(calculate_sds_draws_matrix(dm), expected_res)
})

test_that("statistics_from_single_fit", {
    params <- posterior::as_draws_matrix(
        posterior::as_draws_rvars(list(
            mu = 4,
            tau = 4,
            theta = seq(3.5, 6.5, length.out = 8))))

    # Can't really check correctness, only
    # testing that no error is thrown and structure is OK
    res <- statistics_from_single_fit(posterior::example_draws(example = "eight_schools"),
                               parameters = params, thin_ranks = 1)

    expect_equal(length(unique(res$max_rank)), 1)
    expect_true(all(res$rank >= 0 & res$rank < res$max_rank))
    expect_equal(res$simulated_value, as.numeric(params))
    expect_identical(res$mean > res$simulated_value, sign(res$z_score) < 0)

})
