context("The function that calculates the residuals")

# In this file we test the resid_ID function

test_that("The function gives reasonable results when used correctly by the user", {
  expect_equal(length(resid_ID(rnorm(2000) + c(rep(4, 1000), rep(0, 1000)), cpt = 1000, type_chg = "mean", type_res = "raw")), 2000)
  expect_equal(length(resid_ID(rnorm(1000), cpt = 0, type_chg = "mean", type_res = "standardised")), 1000)
})

test_that("The function gives error messages whenever the input argument is meaningless", {
  # The input for x is a character string.
  expect_error(resid_ID(x = "I am not a numeric vector"))

  # The input for x contains NA
  expect_error(resid_ID(x = c(rnorm(100), NA)))
})
