context("Applying the Isolate-Detect methodology for the case of piecewise-linear mean signals and heavy-tailed noise")

# In this file we test the ht_ID_cplm function

test_that("The correct change-points are given when the function is used correctly", {
  set.seed(111)
  no.cpt <- rt(800, df = 5)
  single.cpt <- c(seq(0, 1999, 1), seq(1998, -1, -1)) + rt(4000, df = 5)
  three.cpts <- c(seq(0, 3998, 2), seq(3996, -2, -2), seq(0, 3998, 2), seq(3996, -2, -2)) + rt(8000, df = 5)

  expect_equal(ht_ID_cplm(no.cpt)$cpt, 0)
  expect_equal(ht_ID_cplm(no.cpt)$no_cpt, 0)

  expect_lt(ht_ID_cplm(single.cpt)$cpt, 2010)
  expect_gt(ht_ID_cplm(single.cpt)$cpt, 1990)
  expect_equal(ht_ID_cplm(single.cpt)$no_cpt, 1)

  expect_lt(ht_ID_cplm(three.cpts)$cpt[1], 2010)
  expect_gt(ht_ID_cplm(three.cpts)$cpt[1], 1990)
  expect_lt(ht_ID_cplm(three.cpts)$cpt[2], 4010)
  expect_gt(ht_ID_cplm(three.cpts)$cpt[2], 3990)
  expect_lt(ht_ID_cplm(three.cpts)$cpt[3], 6010)
  expect_gt(ht_ID_cplm(three.cpts)$cpt[3], 5990)
  expect_equal(ht_ID_cplm(three.cpts)$no_cpt, 3)
})

test_that("Error and warning messages are given correctly", {
  set.seed(111)
  test_data <- rt(1000, df = 5)

  # A string vector is given as the input for xd.
  expect_error(ht_ID_cplm(x = "I am not a numeric"))

  # The threshold constant is set equal to zero.
  expect_error(ht_ID_cplm(x = test_data, ht_thr_id  = 0))
  expect_error(ht_ID_cplm(x = test_data, ht_th_ic_id = -10))

  # The value for lambda is negative.
  expect_error(ht_ID_cplm(x = test_data, p_thr = -2))
  expect_error(ht_ID_cplm(x = test_data, p_ic = -2))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(ht_ID_cplm(x = test_data, p_thr = 3.4))
  expect_warning(ht_ID_cplm(x = test_data, p_ic = 5.4))
})
