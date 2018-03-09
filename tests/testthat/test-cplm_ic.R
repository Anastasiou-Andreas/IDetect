context("The change-points for piecewise-linear mean signals through the information criterion")

# In this file we test the cplm_ic function

test_that("The correct change-points are given when the function is used correctly", {
  set.seed(100)
  no.cpt <- rnorm(1000) + seq(0, 999, 1)
  single.cpt <- rnorm(2000) + c(seq(0, 999, 1), seq(998.5, 499, -0.5))
  three.cpts <- rnorm(2000) + c(seq(0, 499, 1), seq(498.5, 249, -0.5), seq(251, 1249, 2), seq(1248, 749, -1))

  # When there are no change-points, we might get an error message to choose a smaller
  # threshold constant in order to get come estimations that will be used in the solution
  # path.
  expect_equal(cplm_ic(no.cpt, th_const = 1)$cpt_ic$ssic_pen, NA)
  expect_equal(cplm_ic(single.cpt)$cpt_ic$ssic_pen, 1000)
  expect_equal(cplm_ic(three.cpts)$cpt_ic$ssic_pen, c(500, 1000, 1500))
})

test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- rnorm(1000)

  # A string vector is given as the input for xd.
  expect_error(cplm_ic(x = "I am not a numeric"))

  # The threshold constant is set equal to zero.
  expect_error(cplm_ic(x = test_data, th_const = 0))

  # The value for lambda is negative.
  expect_error(cplm_ic(x = test_data, points = -2))

  # The length of x should be greater than 3 in order to have sufficient data
  expect_error(cplm_ic(x = rnorm(3)))

  # Kmax should be positive integer.
  expect_error(cplm_ic(x = test_data, Kmax = 0))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(cplm_ic(x = test_data, points = 3.4))
})
