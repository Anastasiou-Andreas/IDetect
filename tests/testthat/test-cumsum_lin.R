context("The contrast function for the case of piecewise-linear mean signals")

# In this file we test the cumsum_lin function

test_that("The function gives reasonable results when used correctly by the user", {
  expect_equal(length(cumsum_lin(rnorm(100))), 99)
  expect_equal(length(cumsum_lin(rnorm(500))), 499)
  expect_equal(cumsum_lin(rnorm(1)), 0)
  expect_equal(cumsum_lin(rnorm(2)), 0)
})

test_that("The function gives error messages whenever the input argument is meaningless", {
  # The input for x is a character string.
  expect_error(cumsum_lin(x = "I am not a numeric vector"))
})
