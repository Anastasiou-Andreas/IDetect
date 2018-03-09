context("The function that calculates the CUSUM at specific values")

# In this file we test the cusum_one function

test_that("The function gives reasonable results when used correctly by the user", {
  expect_equal(length(cusum_one(rnorm(100), s = c(1, 3, 7, 8), e = c(81, 93, 97, 98), b = c(10, 31, 26, 76))), 4)
  expect_equal(length(cusum_one(rnorm(100), s = c(1, 5), e = c(93, 97), b = c(10, 31))), 2)
  expect_equal(length(cusum_one(rnorm(100), s = 1, e = 100, b = 50)), 1)
})

test_that("The function gives error messages whenever the input argument is meaningless", {
  # The input for x is a character string.
  expect_error(cusum_one(x = "I am not a numeric vector", s = 1, e = 5, b = 3))

  # The vectors below for s, b, e are not of the same length
  expect_error(cusum_one(x = rnorm(100), s = 1, e = c(5, 100), b = c(3, 71)))

  # At least one value of s or b or e is not a positive integer
  expect_error(cusum_one(x = rnorm(100), s = c(-4, 65), e = c(50, 100), b = c(30, 71)))
  expect_error(cusum_one(x = rnorm(100), s = c(4.23, 65), e = c(50, 100), b = c(30, 71)))

  # The value for b[j] should be in the interval [s[j], e[j])
  expect_error(cusum_one(x = rnorm(100), s = c(4, 65), e = c(5, 100), b = c(3, 71)))
})
