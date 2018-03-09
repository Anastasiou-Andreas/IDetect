context("Applying the Isolate-Detect methodology for the case of piecewise-linear mean signals")

# In this file we test the ID_cplm function

test_that("The correct change-points are given when the function is used correctly", {
  set.seed(100)
  no.cpt <- rnorm(1000) + seq(0, 999, 1)
  single.cpt <- rnorm(2000) + c(seq(0, 999, 1), seq(998.5, 499, -0.5))
  three.cpts <- rnorm(2000) + c(seq(0, 499, 1), seq(498.5, 249, -0.5), seq(251, 1249, 2), seq(1248, 749, -1))

  # When there are no change-points, we might get an error message to choose a smaller
  # threshold constant in order to get come estimations that will be used in the solution
  # path.
  expect_equal(ID_cplm(no.cpt, th_ic_id = 0.6)$cpt, 0)
  expect_equal(ID_cplm(no.cpt, th_ic_id = 0.6)$no_cpt, 0)

  expect_equal(ID_cplm(single.cpt)$cpt, 1000)
  expect_equal(ID_cplm(single.cpt)$no_cpt, 1)

  expect_equal(ID_cplm(three.cpts)$cpt, c(500, 1000, 1500))
  expect_equal(ID_cplm(three.cpts)$no_cpt, 3)
})

test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- rnorm(1000)

  # A string vector is given as the input for xd.
  expect_error(ID_cplm(x = "I am not a numeric"))

  # The threshold constant is set equal to zero.
  expect_error(ID_cplm(x = test_data, thr_id  = 0))
  expect_error(ID_cplm(x = test_data, th_ic_id = -10))

  # The value for lambda is negative.
  expect_error(ID_cplm(x = test_data, pointsth = -2))
  expect_error(ID_cplm(x = test_data, pointsic = -2))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(ID_cplm(x = test_data, pointsth = 3.4))
  expect_warning(ID_cplm(x = test_data, pointsic = 5.4))
})
