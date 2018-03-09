context("Applying the Isolate-Detect methodology for the case of piecewise-constant mean signals")

# In this file we test the ID_pcm function

test_that("The correct change-points are given when the function is used correctly", {
  set.seed(100)
  no.cpt <- rnorm(100)
  one.cpt <- rnorm(100) + c(rep(-10, 50), rep(20, 50))
  two.cpts <- rnorm(500) + c(rep(0, 200), rep(10, 150), rep(-10, 150))
  many.cpts <- rnorm(1000) + rep(c(rep(0, 100), rep(10, 100)), 5)

  # When there are no change-points, we will get an error message to choose a smaller
  # threshold constant in order to get come estimations that will be used in the solution
  # path.
  expect_equal(ID_pcm(no.cpt, th_ic_id = 0.6)$cpt, 0)
  expect_equal(ID_pcm(no.cpt, th_ic_id = 0.6)$no_cpt, 0)

  expect_equal(ID_pcm(one.cpt)$cpt, 50)
  expect_equal(ID_pcm(one.cpt)$no_cpt, 1)

  expect_equal(ID_pcm(two.cpts)$cpt, c(200, 350))
  expect_equal(ID_pcm(two.cpts)$no_cpt, 2)

  expect_equal(ID_pcm(many.cpts)$cpt, seq(100, 900, 100))
  expect_equal(ID_pcm(many.cpts)$no_cpt, 9)
})

test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- rnorm(1000)

  # A string vector is given as the input for xd.
  expect_error(ID_pcm(x = "I am not a numeric"))

  # The threshold constant is set equal to zero.
  expect_error(ID_pcm(x = test_data, thr_id  = 0))
  expect_error(ID_pcm(x = test_data, th_ic_id = -10))

  # The value for lambda is negative.
  expect_error(ID_pcm(x = test_data, pointsth = -2))
  expect_error(ID_pcm(x = test_data, pointsic = -2))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(ID_pcm(x = test_data, pointsth = 3.4))
  expect_warning(ID_pcm(x = test_data, pointsic = 5.4))
})
