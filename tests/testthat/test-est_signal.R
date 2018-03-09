context("The function that estimates the underlying signal of a given data sequence")

# In this file we test the est_signal function

test_that("The function gives reasonable results when used correctly by the user", {
  #no change-points
  expect_equal(length(est_signal(rnorm(100), cpt = 0, type = "mean")), 100)
  expect_equal(length(est_signal(rnorm(500), cpt = 0, type = "mean")), 500)

  #one change-point at position 250
  expect_equal(length(est_signal(rnorm(500) + c(rep(0, 250), rep(10, 250)), cpt = 250, type = "mean")), 500)

  #one change-point at the position 1000 in a piecewise-linear mean signal
  expect_equal(length(est_signal(rnorm(2000) + c(seq(0, 999, 1), seq(998.5, 499, -0.5)), cpt = 1000, type = "slope")), 2000)
})

test_that("The function gives error messages whenever the input argument is meaningless", {
  # The input for x is a character string.
  expect_error(est_signal(x = "I am not a numeric vector"))

  # The change-points cannot be negative or greater than the length of the data sequence.
  expect_error(est_signal(x = rnorm(100), cpt = c(-1, 60), type = "mean"))
  expect_error(est_signal(x = rnorm(100), cpt = c(4, 103), type = "mean"))

  # The input for x cannot have NAs.
  expect_error(est_signal(x = c(rnorm(100), NA), cpt = 50, type = "mean"))
})
