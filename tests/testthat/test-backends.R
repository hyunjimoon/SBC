test_that("draws_rvars_to_standata converts properly", {

    scalar_var <- 1
    vector_var <- array(1:10, 10)
    matrix_var <- array(1:12, c(3,4))

    list_val <- list(scalar_var = scalar_var,
                     vector_var = vector_var,
                     matrix_var = matrix_var)

    rvar_val <- posterior::draws_rvars(scalar_var = posterior::rvar(scalar_var),
                            vector_var = posterior::rvar(array(vector_var, dim = c(1, dim(vector_var))), dim = dim(vector_var)),
                            matrix_var = posterior::rvar(array(matrix_var, dim = c(1, dim(matrix_var))), dim = dim(matrix_var)))

    expect_equal(draws_rvars_to_standata_single(rvar_val), list_val)
})

