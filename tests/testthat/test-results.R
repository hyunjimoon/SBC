test_that("capture_all_outputs", {
    expect_identical(
        capture_all_outputs({
            cat("Test\n")
            warning("W")
            message("M", appendLF = FALSE)
            warning("W2")
            message("M2", appendLF = FALSE)
            message("M3", appendLF = FALSE)

            # A special case - silent error
            try(stop("Error"))

            14
        }),
        list(result = 14,
             messages = c("M", "M2", "M3"),
             warnings = c("W", "W2"),
             output = c('Test', 'Error in try(stop("Error")) : Error')))

    # Nested capture.output

    expect_identical(
        capture_all_outputs({
            captured <- capture_all_outputs({
                cat("Test\n")
                warning("W")
                message("M", appendLF = FALSE)

                # A special case - silent error
                try(stop("Error"))

                28

            })
            cat("BEFORE\n")
            message("M_BEFORE", appendLF = FALSE)
            warning("W_BEFORE")
            try(stop("E_BEFORE"))
            reemit_captured(captured)
            try(stop("E_AFTER"))
            warning("W_AFTER")
            message("M_AFTER", appendLF = FALSE)
            cat("AFTER\n")
            13
        }),
        list(result = 13,
             messages = c("M_BEFORE", "M", "M_AFTER"),
             warnings = c("W_BEFORE", "W", "W_AFTER"),
             output = c('BEFORE',
                        'Error in try(stop("E_BEFORE")) : E_BEFORE',
                        'Test',
                        'Error in try(stop("Error")) : Error',
                        'Error in try(stop("E_AFTER")) : E_AFTER',
                        'AFTER'
                        ))

    )
})

test_that("subset_bind", {
    res <- SBC_results(stats = data.frame(sim_id = rep(1:3, each = 4), s = 1:12),
                       fits = list("A", NULL, "C"),
                       outputs = list(c("A1","A2"), c(), c("C1", "C4")),
                       warnings = list(c(), "XXXX", "asdfdaf"),
                       messages = list("aaaa", "ddddd", NA_character_),
                       errors = list(NULL, "customerror", NULL),
                       default_diagnostics = data.frame(sim_id = 1:3, qq = rnorm(3)),
                       backend_diagnostics = data.frame(sim_id = 1:3, rr = rnorm(3))
                       )

    remove_sim_id_names <- function(x) {
        names(x$stats$sim_id) <- NULL
        names(x$default_diagnostics$sim_id) <- NULL
        if(!is.null(x$backend_diagnostics)) {
          names(x$backend_diagnostics$sim_id) <- NULL
        }
        x
    }

    expect_equal(res, remove_sim_id_names(bind_results(res[1], res[2:3])))
    expect_equal(res, remove_sim_id_names(bind_results(res[1:2], res[3])))
    expect_equal(remove_sim_id_names(res[3:1]), remove_sim_id_names(bind_results(res[3:2], res[1])))
    expect_equal(remove_sim_id_names(res[2]), remove_sim_id_names(((res[2:3])[1])))

    # The same, but with some NULLs
    res2 <- SBC_results(stats = data.frame(sim_id = rep(1:3, each = 4), s = 1:12),
                       fits = list("A", NULL, "C"),
                       outputs = rep(list(NULL), 3),
                       warnings = rep(list(NULL), 3),
                       messages = rep(list(NULL), 3),
                       errors = rep(list(NULL), 3),
                       default_diagnostics = data.frame(sim_id = 1:3, qq = rnorm(3)),
                       backend_diagnostics = NULL
    )

    expect_equal(res2, remove_sim_id_names(bind_results(res2[1], res2[2:3])))
    expect_equal(res2, remove_sim_id_names(bind_results(res2[1:2], res2[3])))
    expect_equal(remove_sim_id_names(res2[3:1]), remove_sim_id_names(bind_results(res2[3:2], res2[1])))
    expect_equal(remove_sim_id_names(res2[2]), remove_sim_id_names(((res2[2:3])[1])))


})

test_that("calculate_ranks_draws_matrix works", {

    dm <- matrix(NA_real_, nrow = 10, ncol = 4)
    colnames(dm) <- c("a","b","c", "d")

    dm[, "a"] <- sample(1:10)
    dm[, "b"] <- sample(1:10)
    dm[, "c"] <- sample(1:10)
    dm[, "d"] <- sample(c(1:5, 1:5))
    dm <- posterior::as_draws_matrix(dm)

    vars <- matrix(c(3.5, -5, 15, 3), nrow = 1)
    colnames(vars) <- c("a","b","c", "d")

    N_steps <- 1e4
    all_ranks <- matrix(NA_real_, nrow = N_steps, ncol = 4)
    for(i in 1:N_steps) {
        last_ranks <- calculate_ranks_draws_matrix(vars, dm)
        all_ranks[i,] <- last_ranks

    }
    expect_true(!any(is.na(all_ranks)))
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

test_that("calculate_ranks_draws_matrix infinity NA", {
  dm <- matrix(NA_real_, nrow = 5, ncol = 6)
  colnames(dm) <- c("a","b", "c","d", "e", "f")
  dm[,"a"] <- c(-Inf, -Inf, 1, 2, +Inf)
  dm[,"b"] <- c(-Inf, -Inf, -Inf, 1, 2)
  dm[,"c"] <- c(-Inf, NA_real_, -Inf, 1, 2)
  dm[,"d"] <- rep(NA_real_, 5)
  dm[,"e"] <- rep(-Inf, 5)
  dm[,"f"] <- c(-Inf, NA_real_, Inf, -300, -500)

  dm <- posterior::as_draws_matrix(dm)

  vars <- matrix(c(-Inf, NA_real_, 14, NA_real_, -Inf, NA_real_), nrow = 1)
  colnames(vars) <- c("a","b","c", "d", "e", "f")

  # First with default settings
  N_steps <- 200
  all_ranks <- matrix(NA_real_, nrow = N_steps, ncol = ncol(dm))
  for(i in 1:N_steps) {
    last_ranks <- calculate_ranks_draws_matrix(vars, dm)
    all_ranks[i,] <- last_ranks

  }
  expect_true(!any(is.na(all_ranks)))

  expect_true(all(all_ranks[,1] <= 2))
  expect_true(all(0:2 %in% all_ranks[,1]))

  expect_true(all(0:5 %in% all_ranks[,2]))
  expect_true(all(all_ranks[,3] <= 5 & all_ranks[,3] >= 4))
  expect_true(all(4:5 %in% all_ranks[,3]))
  expect_true(all(0:5 %in% all_ranks[,4]))
  expect_true(all(0:5 %in% all_ranks[,5]))
  expect_true(all(0:1 %in% all_ranks[,6]))

  # Now with na_lowest = TRUE
  N_steps <- 200
  all_ranks <- matrix(NA_real_, nrow = N_steps, ncol = ncol(dm))
  for(i in 1:N_steps) {
    last_ranks <- calculate_ranks_draws_matrix(vars, dm, na_lowest = TRUE)
    all_ranks[i,] <- last_ranks

  }
  expect_true(!any(is.na(all_ranks)))

  expect_true(all(all_ranks[,1] <= 2))
  expect_true(all(0:2 %in% all_ranks[,1]))
  expect_true(all(all_ranks[,2] == 0))
  expect_true(all(all_ranks[,3] == 5))
  expect_true(all(0:5 %in% all_ranks[,4]))
  expect_true(all(0:5 %in% all_ranks[,5]))
  expect_true(all(0:1 %in% all_ranks[,6]))
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

test_that("SBC_statistics_from_single_fit", {
    vars <- posterior::as_draws_matrix(
        posterior::draws_rvars(
            mu = posterior::rvar(4) ,
            tau = posterior::rvar(4),
            theta = posterior::rvar(array(seq(3.5, 6.5, length.out = 8), dim = c(1,8)))))

    # Can't really check correctness, only
    # testing that no error is thrown and structure is OK
    test_draws <- posterior::example_draws(example = "eight_schools")
    res <- SBC_statistics_from_single_fit(test_draws,
                               variables = vars, thin_ranks = 1, dquants = NULL,
                               ensure_num_ranks_divisor = 1,
                               backend = SBC_backend_mock(),
                               var_attributes = NULL)


    expect_equal(length(unique(res$max_rank)), 1)
    expect_equal(unique(res$max_rank), posterior::ndraws(test_draws))
    expect_true(all(res$rank >= 0 & res$rank < res$max_rank))
    expect_equal(res$simulated_value, as.numeric(vars))
    expect_identical(res$mean > res$simulated_value, sign(res$z_score) < 0)

    # Test ensure_num_ranks_divisor
    # Make sure the test draws have the expected size before proceeding
    expect_equal(posterior::ndraws(test_draws), 400)
    res_ensure2 <- SBC_statistics_from_single_fit(posterior::example_draws(example = "eight_schools"),
                                      variables = vars, thin_ranks = 1, dquants = NULL,
                                      ensure_num_ranks_divisor = 2,
                                      backend = SBC_backend_mock(),
                                      var_attributes = NULL)
    # Number of ranks = max_rank + 1 (as 0 is a valid rank)
    expect_equal(unique(res_ensure2$max_rank), 399)


    # Test ensure_num_ranks_divisor, combined with thin_ranks
    res_ensure7 <- SBC_statistics_from_single_fit(posterior::example_draws(example = "eight_schools"),
                                              variables = vars, thin_ranks = 4, dquants = NULL,
                                              ensure_num_ranks_divisor = 7,
                                              backend = SBC_backend_mock(),
                                              var_attributes = NULL)
    expect_equal(unique(res_ensure7$max_rank), 97)

})

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
