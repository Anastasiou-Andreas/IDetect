context("The change-points for piecewise-constant mean through the information criterion")

# In this file we test the pcm_ic function

test_that("The correct change-points are given when the function is used correctly", {
  set.seed(100)
  no.cpt <- rnorm(100)
  one.cpt <- rnorm(100) + c(rep(-10, 50), rep(20, 50))
  two.cpts <- rnorm(500) + c(rep(0, 200), rep(10, 150), rep(-10, 150))
  many.cpts <- rnorm(1000) + rep(c(rep(0, 100), rep(10, 100)), 5)

  # When there are no change-points, we will get an error message to choose a smaller
  # threshold constant in order to get come estimations that will be used in the solution
  # path.
  expect_equal(pcm_ic(no.cpt, th_const = 0.6)$cpt_ic$ssic_pen, NA)
  expect_equal(pcm_ic(one.cpt)$cpt_ic$ssic_pen, 50)
  expect_equal(pcm_ic(two.cpts)$cpt_ic$ssic_pen, c(200, 350))
  expect_equal(pcm_ic(many.cpts)$cpt_ic$ssic_pen, seq(100, 900, 100))
            })

test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- rnorm(1000)

  # A string vector is given as the input for xd.
  expect_error(pcm_ic(x = "I am not a numeric"))

  # The threshold constant is set equal to zero.
  expect_error(pcm_ic(x = test_data, th_const = 0))

  # The value for lambda is negative.
  expect_error(pcm_ic(x = test_data, points = -2))

  # The length of x should be greater than 3 in order to have sufficient data
  expect_error(pcm_ic(x = rnorm(3)))

  # Kmax should be positive integer.
  expect_error(pcm_ic(x = test_data, Kmax = 0))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(pcm_ic(x = test_data, points = 3.4))
})
