context("The function that calculates the CUSUM")

# In this file we test the function with the name cusum_function

test_that("The function gives reasonable results when used correctly by the user", {
  expect_equal(length(cusum_function(rnorm(100))), 99)
  expect_equal(length(cusum_function(rnorm(500))), 499)

  # When x has length equal to one, then the result will be the vector (NaN, -Inf)
  only.one.pt <- cusum_function(rnorm(1))
  expect_equal(length(only.one.pt), 2)
  expect_equal(only.one.pt[1], NaN)
  expect_equal(abs(only.one.pt[2]), Inf)
})

test_that("The function gives error messages whenever the input argument is meaningless", {
  # The input for x is a character string.
  expect_error(cusum_function(x = "I am not a numeric vector"))
})
