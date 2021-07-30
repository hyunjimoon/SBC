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

