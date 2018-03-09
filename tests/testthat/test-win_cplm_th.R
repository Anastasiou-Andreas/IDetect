context("The change-points for piecewise-linear mean through thresholding and using the windows structure")

# In this file we test the win_cplm_th function

test_that("The correct change-points are given through the function cplm_th
          as the length of x is relatively small and the windows structure is not needed", {
  set.seed(100)
  no.cpt <- rnorm(1000) + seq(0, 999, 1)
  single.cpt <- rnorm(2000) + c(seq(0, 999, 1), seq(998.5, 499, -0.5))
  three.cpts <- rnorm(2000) + c(seq(0, 499, 1), seq(498.5, 249, -0.5), seq(251, 1249, 2), seq(1248, 749, -1))
  expect_equal(win_cplm_th(no.cpt), numeric(0))
  expect_equal(win_cplm_th(single.cpt), 1000)
  expect_equal(win_cplm_th(three.cpts), c(500, 1000, 1500))
})

test_that("The correct change-points are given through the windows approach, which is mainly employed for
          larger time series", {
  set.seed(101)
  win.three.cpts <- rnorm(25000) + c(seq(0, 5999, 1), seq(5998.5, 2999, -0.5), seq(3001, 14999, 2), seq(14998, 7999, -1))
  win.no.cpt <- rnorm(15000)
  expect_equal(win_cplm_th(win.three.cpts), c(6000, 12000, 18000))
  expect_equal(win_cplm_th(win.no.cpt), numeric(0))
            })

test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- rnorm(15000)

  # A string vector is given as the input for xd.
  expect_error(win_cplm_th(xd = "I am not a numeric"))

  # The threshold constant is set equal to zero.
  expect_error(win_cplm_th(xd = test_data, thr_con = 0))

  # The value for lambda is negative.
  expect_error(win_cplm_th(xd = test_data, w_points = -2))

  # The value for the windows length is negative.
  expect_error(win_cplm_th(xd = test_data, c_win = -2))

  # The value for l_win is negative.
  expect_error(win_cplm_th(xd = test_data, l_win = -200))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(win_cplm_th(xd = test_data, w_points = 3.4))

  # The value for c_win is a positive real number instead of a positive integer.
  expect_warning(win_cplm_th(xd = test_data, c_win = 103.4))

  # The value for l_win is a positive real number instead of a positive integer.
  expect_warning(win_cplm_th(xd = test_data, l_win = 2000.4))
})
