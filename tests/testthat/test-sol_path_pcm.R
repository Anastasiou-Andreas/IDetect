context("The function that calculates the solution path in the case of piecewise-constant mean")

# In this file we test the sol_path_pcm function

test_that("To test that the solution path includes the true change-points which should appear before
          any other values for the change-points", {
  set.seed(100)
  one.cpt <- rnorm(100) + c(rep(-10, 50), rep(20, 50))
  two.cpts <- rnorm(500) + c(rep(0, 200), rep(10, 150), rep(-10, 150))
  three.cpts <- rnorm(1000) + c(rep(0, 250), rep(10, 250), rep(-10, 250), rep(10, 250))
  sol.p1 <- sol_path_pcm(one.cpt)
  expect_equal(sol.p1[1], 50)

  sol.p2 <- sol_path_pcm(two.cpts)
  expect_true( (sol.p2[1] == 200) || (sol.p2[1] == 350))
  expect_true( ( (sol.p2[1] == 200) || (sol.p2[1] == 350)) && (sol.p2[1] != sol.p2[2]))

  sol.p3 <- sol_path_pcm(three.cpts)
  expect_true( (sol.p3[1] == 250) || (sol.p3[1] == 500) || (sol.p3[1] == 750))
  expect_true( ( (sol.p3[2] == 250) || (sol.p3[2] == 500) || (sol.p3[2] == 750)) && (sol.p3[1] != sol.p3[2]))
  expect_true( ( (sol.p3[3] == 250) || (sol.p3[3] == 500) || (sol.p3[3] == 750)) && (sol.p3[3] != sol.p3[2]) && (sol.p3[3] != sol.p3[1]))
})

test_that("The function gives error messages whenever the input argument is meaningless", {
  # The input for x is a character string.
  expect_error(sol_path_pcm(x = "I am not a numeric vector"))

  # The threshold constant is set equal to zero.
  expect_error(sol_path_pcm(x = rnorm(100), thr_ic = 0))

  # The value for lambda is negative.
  expect_error(sol_path_pcm(x = rnorm(100), points = -2))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(sol_path_pcm(x = rnorm(100), points = 3.4))
})
