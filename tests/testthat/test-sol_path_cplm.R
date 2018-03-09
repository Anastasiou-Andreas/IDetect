context("The function that calculates the solution path in the case of piecewise-linear mean")

# In this file we test the sol_path_cplm function

test_that("To test that the solution path includes the true change-points which should appear before
          any other values for the change-points", {
  set.seed(300)
  single.cpt <- rnorm(2000) + c(seq(0, 999, 1), seq(998.5, 499, -0.5))
  two.cpts <- rnorm(2000) + c(seq(0, 699, 1), seq(698.5, 399, -0.5), seq(401, 1799, 2))
  three.cpts <- rnorm(2000) + c(seq(0, 499, 1), seq(498.5, 249, -0.5), seq(251, 1249, 2), seq(1248, 749, -1))
  sol.p1 <- sol_path_cplm(single.cpt)
  expect_equal(sol.p1[1], 1000)

  sol.p2 <- sol_path_cplm(two.cpts)
  expect_true( (sol.p2[1] == 700) || (sol.p2[1] == 1300))
  expect_true( ( (sol.p2[1] == 700) || (sol.p2[1] == 1300)) && (sol.p2[1] != sol.p2[2]))
})

test_that("The function gives error messages whenever the input argument is meaningless", {
  # The input for x is a character string.
  expect_error(sol_path_cplm(x = "I am not a numeric vector"))

  # The threshold constant is set equal to zero.
  expect_error(sol_path_cplm(x = rnorm(100), thr_ic = 0))

  # The value for lambda is negative.
  expect_error(sol_path_cplm(x = rnorm(100), points = -2))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(sol_path_cplm(x = rnorm(100), points = 3.4))
})
