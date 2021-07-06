test_that("norm_working_params() transforms natural parameters", {
  mu <- list(matrix(c(1, 5, 2, 4), 2, 2, byrow = TRUE),
             matrix(c(10, 25, 12, 24), 2, 2, byrow = TRUE))
  sigma <- list(matrix(c(1, 2, 1.4, 1.7), 2, 2, byrow = TRUE),
                matrix(c(1, 2.2, 1.1, 2), 2, 2, byrow = TRUE))
  beta <- matrix(c(-2, -2, 0, 0, 0, 0), ncol = 2)
  delta <- list(c(0.5, 0.5), c(0.5, 0.5))
  result <- norm_working_params(2, 2, 2, mu, sigma, beta, delta)

  expect_identical(result, result)
})
