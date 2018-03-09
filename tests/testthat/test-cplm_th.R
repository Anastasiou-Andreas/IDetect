context("The change-points for piecewise-linear mean through thresholding")

# In this file we test the cplm_th function

test_that("The correct change-points are given", {
  set.seed(100)
  no.cpt <- rnorm(1000) + seq(0, 999, 1)
  single.cpt <- rnorm(2000) + c(seq(0, 999, 1), seq(998.5, 499, -0.5))
  three.cpts <- rnorm(2000) + c(seq(0, 499, 1), seq(498.5, 249, -0.5), seq(251, 1249, 2), seq(1248, 749, -1))
  expect_equal(cplm_th(no.cpt), numeric(0))
  expect_equal(cplm_th(single.cpt), 1000)
  expect_equal(cplm_th(three.cpts), c(500, 1000, 1500))
})

test_that("Error and warning messages are given correctly", {
  # A string vector is given as the input for x.
  expect_error(cplm_th(x = "I am not a numeric"))

  # The threshold constant is equal to zero.
  expect_error(cplm_th(x = rnorm(1000), thr_const = 0))

  # The value for lambda is negative.
  expect_error(cplm_th(x = rnorm(1000), points = -2))

  # The value for lambda is a positive real number, instead of a positive integer
  expect_warning(cplm_th(x = rnorm(1000), points = 3.4))
})
