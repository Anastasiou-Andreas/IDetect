context("Applying the Isolate-Detect methodology for the case of piecewise-constant mean signals with heavy-tailed noise")

# In this file we test the ht_ID_pcm function

test_that("The correct change-points are given when the function is used correctly", {
  set.seed(10)
  no.cpt <- rt(1000, df = 5)
  single.cpt <- c(rep(4, 3000), rep(0, 3000)) + rt(6000, df = 5)
  three.cpts <- c(rep(4, 2000), rep(0, 2000), rep(-4, 2000), rep(0, 2000)) + rt(8000, df = 5)

  # When there are no change-points, we will get an error message to choose a smaller
  # threshold constant in order to get come estimations that will be used in the solution
  # path.
  expect_equal(ht_ID_pcm(no.cpt, ht_th_ic_id = 0.9)$cpt, 0)
  expect_equal(ht_ID_pcm(no.cpt, ht_th_ic_id = 0.9)$no_cpt, 0)

  expect_lt(ht_ID_pcm(single.cpt)$cpt, 3010)
  expect_gt(ht_ID_pcm(single.cpt)$cpt, 2990)
  expect_equal(ht_ID_pcm(single.cpt)$no_cpt, 1)

  expect_lt(ht_ID_pcm(three.cpts)$cpt[1], 2010)
  expect_gt(ht_ID_pcm(three.cpts)$cpt[1], 1990)
  expect_lt(ht_ID_pcm(three.cpts)$cpt[2], 4010)
  expect_gt(ht_ID_pcm(three.cpts)$cpt[2], 3990)
  expect_lt(ht_ID_pcm(three.cpts)$cpt[3], 6010)
  expect_gt(ht_ID_pcm(three.cpts)$cpt[3], 5990)
  expect_equal(ht_ID_pcm(three.cpts)$no_cpt, 3)
})

test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- rt(1000, df = 5)

  # A string vector is given as the input for xd.
  expect_error(ht_ID_pcm(x = "I am not a numeric"))

  # The threshold constant is set equal to zero.
  expect_error(ht_ID_pcm(x = test_data, ht_thr_id  = 0))
  expect_error(ht_ID_pcm(x = test_data, ht_th_ic_id = -10))

  # The value for lambda is negative.
  expect_error(ht_ID_pcm(x = test_data, p_thr = -2))
  expect_error(ht_ID_pcm(x = test_data, p_ic = -2))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(ht_ID_pcm(x = test_data, p_thr = 3.4, ht_th_ic_id = 0.6))
  expect_warning(ht_ID_pcm(x = test_data, p_ic = 5.4, ht_th_ic_id = 0.6))
})
