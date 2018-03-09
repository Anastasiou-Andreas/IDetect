context("The function that transforms the given sequence in order to have noise closer to the Gaussian distribution")

# In this file we test the normalise function

test_that("The obtained sequence has the correct length", {
  expect_equal(length(normalise(rt(n = 10000, df = 5), sc = 3)), 3334)

  # When the value for sc is greater that the length of the given data sequence x
  # then the result is the mean of x times the square root of n, where n is the length
  # of x
  t1 <- rt(n = 10, df = 5)
  expect_equal(normalise(t1, sc = 100), sqrt(10) * mean(t1))
  expect_equal(length(normalise(t1, sc = 100)), 1)
})

test_that("The function gives error messages whenever the input argument is meaningless", {
  # The input for x is a character string.
  expect_error(normalise(x = "I am not a numeric vector"))

  # The input for scale should be a positive integer.
  expect_error(normalise(x = rt(100, df = 5), sc = -3))

  # The value for sc is a positive real number instead of a positive integer.
  expect_warning(normalise(x = rt(100, df = 5), sc = 3.8))
})
