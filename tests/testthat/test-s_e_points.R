context("Start- and end-points")

# In this file we test the s_e_points function

test_that("s_e_points works for reasonable inputs by the user", {
  example1 <- s_e_points(r = seq(10, 100, 10), l = seq(91, 1, -10), s = 5, e = 32)
  expect_equal(example1[[1]], c(10, 20, 30, 32))
  expect_equal(example1[[2]], c(31, 21, 11, 5))

  example2 <- s_e_points(r = seq(5, 200, 5), l = seq(196, 1, -5), s = 87, e = 112)
  expect_equal(example2[[1]], c(90, 95, 100, 105, 110, 112))
  expect_equal(example2[[2]], c(111, 106, 101, 96, 91, 87))

  example3 <- s_e_points(r = seq(2, 20, 1), l = seq(19, 1, -1), s = 13, e = 16)
  expect_equal(example3[[1]], c(14, 15, 16))
  expect_equal(example3[[2]], c(15, 14, 13))
})

test_that("s_e_points gives an error message whenever needed", {
  #r and l have negative entries
  expect_error(s_e_points(r = seq(-10, 100, 10), l = seq(99, -9, -10), s = 5, e = 32))

  #s is negative
  expect_error(s_e_points(r = seq(10, 100, 10), l = seq(99, 1, -10), s = -5, e = 32))

  #s is greater than e
  expect_error(s_e_points(r = seq(10, 100, 10), l = seq(99, 1, -10), s = 25, e = 12))
})

test_that("s_e_points gives an error message whenever needed", {
  #s or/and e are positive real numbers and their integer part is used in order to evaluate the function
  expect_warning(s_e_points(r = seq(10, 100, 10), l = seq(91, 1, -10), s = 5.6, e = 32.4))
  example4 <- s_e_points(r = seq(10, 100, 10), l = seq(91, 1, -10), s = 5.6, e = 32.4)
  expect_equal(example4[[1]], c(10, 20, 30, 32))
  expect_equal(example4[[2]], c(31, 21, 11, 5))
})
