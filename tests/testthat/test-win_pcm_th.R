context("The change-points for piecewise-constant mean through thresholding and using the windows structure")

# In this file we test the win_pcm_th function

test_that("The correct change-points are given through the function pcm_th
          as the length of x is relatively small and the windows structure is not needed", {
  set.seed(100)
  no.cpt <- rnorm(100)
  one.cpt <- rnorm(100) + c(rep(-10, 50), rep(20, 50))
  two.cpts <- rnorm(500) + c(rep(0, 200), rep(10, 150), rep(-10, 150))
  many.cpts <- rnorm(1000) + rep(c(rep(0, 100), rep(10, 100)), 5)
  expect_equal(win_pcm_th(no.cpt), numeric(0))
  expect_equal(win_pcm_th(one.cpt), 50)
  expect_equal(win_pcm_th(two.cpts), c(200, 350))
  expect_equal(win_pcm_th(many.cpts), seq(100, 900, 100))
})

test_that("The correct change-points are given through the windows approach, which is mainly employed for
          larger time series", {
  set.seed(101)
  win.many.cpts <- rnorm(22500) + rep(c(rep(0, 1500), rep(10, 3000)), 5)
  win.no.cpt <- rnorm(15000)
  expect_equal(win_pcm_th(win.many.cpts), c(1500, 4500, 6000, 9000, 10500, 13500, 15000, 18000, 19500))
  expect_equal(win_pcm_th(win.no.cpt), numeric(0))
})

test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- rnorm(15000)

  # A string vector is given as the input for xd.
  expect_error(win_pcm_th(xd = "I am not a numeric"))

  # The threshold constant is set equal to zero.
  expect_error(win_pcm_th(xd = test_data, thr_con = 0))

  # The value for lambda is negative.
  expect_error(win_pcm_th(xd = test_data, w_points = -2))

  # The value for the windows length is negative.
  expect_error(win_pcm_th(xd = test_data, c_win = -2))

  # The value for l_win is negative.
  expect_error(win_pcm_th(xd = test_data, l_win = -200))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(win_pcm_th(xd = test_data, w_points = 3.4))

  # The value for c_win is a positive real number instead of a positive integer.
  expect_warning(win_pcm_th(xd = test_data, c_win = 103.4))

  # The value for l_win is a positive real number instead of a positive integer.
  expect_warning(win_pcm_th(xd = test_data, l_win = 2000.4))
})
