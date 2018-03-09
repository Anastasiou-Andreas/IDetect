context("Applying the Isolate-Detect methodology")

# In this file we test the ID function

test_that("The correct change-points are given when the function is used correctly
          in the case of piecewise-constant mean signals with Gaussian noise", {
  set.seed(111)
  no.cpt <- rnorm(100)
  one.cpt <- rnorm(100) + c(rep(-10, 50), rep(20, 50))
  two.cpts <- rnorm(500) + c(rep(0, 200), rep(10, 150), rep(-10, 150))
  many.cpts <- rnorm(1000) + rep(c(rep(0, 100), rep(10, 100)), 5)

  # When there are no change-points, we will get an error message to choose a smaller
  # threshold constant in order to get come estimations that will be used in the solution
  # path.
  expect_equal(ID(no.cpt, th.ic = 0.6, contrast = "mean", ht = FALSE)$cpt, 0)
  expect_equal(ID(no.cpt, th.ic = 0.6, contrast = "mean", ht = FALSE)$no_cpt, 0)

  expect_equal(ID(one.cpt, contrast = "mean", ht = FALSE)$cpt, 50)
  expect_equal(ID(one.cpt, contrast = "mean", ht = FALSE)$no_cpt, 1)

  expect_equal(ID(two.cpts, contrast = "mean", ht = FALSE)$cpt, c(200, 350))
  expect_equal(ID(two.cpts, contrast = "mean", ht = FALSE)$no_cpt, 2)

  expect_equal(ID(many.cpts, contrast = "mean", ht = FALSE)$cpt, seq(100, 900, 100))
  expect_equal(ID(many.cpts, contrast = "mean", ht = FALSE)$no_cpt, 9)
})

test_that("The correct change-points are given when the function is used correctly in the
          case of piecewise-constant signals with heavy-tailed noise", {
  set.seed(10)
  no.cpt <- rt(1000, df = 5)
  single.cpt <- c(rep(4, 3000), rep(0, 3000)) + rt(6000, df = 5)
  three.cpts <- c(rep(4, 2000), rep(0, 2000), rep(-4, 2000), rep(0, 2000)) + rt(8000, df = 5)

  # When there are no change-points, we will get an error message to choose a smaller
  # threshold constant in order to get come estimations that will be used in the solution
  # path.
  expect_equal(ID(no.cpt, th.ic = 0.6, contrast = "mean", ht = TRUE)$cpt, 0)
  expect_equal(ID(no.cpt, th.ic = 0.6, contrast = "mean", ht = TRUE)$no_cpt, 0)

  expect_lt(ID(single.cpt, contrast = "mean", ht = TRUE)$cpt, 3010)
  expect_gt(ID(single.cpt, contrast = "mean", ht = TRUE)$cpt, 2990)
  expect_equal(ID(single.cpt, contrast = "mean", ht = TRUE)$no_cpt, 1)

  expect_lt(ID(three.cpts, contrast = "mean", ht = TRUE)$cpt[1], 2010)
  expect_gt(ID(three.cpts, contrast = "mean", ht = TRUE)$cpt[1], 1990)
  expect_lt(ID(three.cpts, contrast = "mean", ht = TRUE)$cpt[2], 4010)
  expect_gt(ID(three.cpts, contrast = "mean", ht = TRUE)$cpt[2], 3990)
  expect_lt(ID(three.cpts, contrast = "mean", ht = TRUE)$cpt[3], 6010)
  expect_gt(ID(three.cpts, contrast = "mean", ht = TRUE)$cpt[3], 5990)
  expect_equal(ID(three.cpts, contrast = "mean", ht = TRUE)$no_cpt, 3)
})

test_that("The correct change-points are given when the function is used correctly in the
          case of piecewise-linear mean signals with Gaussian noise", {
  set.seed(100)
  no.cpt <- rnorm(1000) + seq(0, 999, 1)
  single.cpt <- rnorm(2000) + c(seq(0, 999, 1), seq(998.5, 499, -0.5))
  three.cpts <- rnorm(2000) + c(seq(0, 499, 1), seq(498.5, 249, -0.5), seq(251, 1249, 2), seq(1248, 749, -1))

  # When there are no change-points, we might get an error message to choose a smaller
  # threshold constant in order to get come estimations that will be used in the solution
  # path.
  expect_equal(ID(no.cpt, th.ic.lin = 0.6, contrast = "slope", ht = FALSE)$cpt, 0)
  expect_equal(ID(no.cpt, th.ic.lin = 0.6, contrast = "slope", ht = FALSE)$no_cpt, 0)

  expect_equal(ID(single.cpt, contrast = "slope", ht = FALSE)$cpt, 1000)
  expect_equal(ID(single.cpt, contrast = "slope", ht = FALSE)$no_cpt, 1)

  expect_equal(ID(three.cpts, contrast = "slope", ht = FALSE)$cpt, c(500, 1000, 1500))
  expect_equal(ID(three.cpts, contrast = "slope", ht = FALSE)$no_cpt, 3)
})

test_that("The correct change-points are given when the function is used correctly in the
          case of piecewise-linear mean with heavy-tailed noise", {
  set.seed(111)
  no.cpt <- rt(800, df = 5)
  single.cpt <- c(seq(0, 1999, 1), seq(1998, -1, -1)) + rt(4000, df = 5)
  three.cpts <- c(seq(0, 3998, 2), seq(3996, -2, -2), seq(0, 3998, 2), seq(3996, -2, -2)) + rt(8000, df = 5)

  expect_equal(ID(no.cpt, contrast = "slope", ht = TRUE)$cpt, 0)
  expect_equal(ID(no.cpt, contrast = "slope", ht = TRUE)$no_cpt, 0)

  expect_lt(ID(single.cpt, contrast = "slope", ht = TRUE)$cpt, 2010)
  expect_gt(ID(single.cpt, contrast = "slope", ht = TRUE)$cpt, 1990)
  expect_equal(ID(single.cpt, contrast = "slope", ht = TRUE)$no_cpt, 1)

  expect_lt(ID(three.cpts, contrast = "slope", ht = TRUE)$cpt[1], 2010)
  expect_gt(ID(three.cpts, contrast = "slope", ht = TRUE)$cpt[1], 1990)
  expect_lt(ID(three.cpts, contrast = "slope", ht = TRUE)$cpt[2], 4010)
  expect_gt(ID(three.cpts, contrast = "slope", ht = TRUE)$cpt[2], 3990)
  expect_lt(ID(three.cpts, contrast = "slope", ht = TRUE)$cpt[3], 6010)
  expect_gt(ID(three.cpts, contrast = "slope", ht = TRUE)$cpt[3], 5990)
  expect_equal(ID(three.cpts, contrast = "slope", ht = TRUE)$no_cpt, 3)
})


test_that("Error and warning messages are given correctly", {
  set.seed(100)
  test_data <- rnorm(1000)

  # A string vector is given as the input for xd.
  expect_error(ID(x = "I am not a numeric"))

  # The threshold constant is set equal to zero.
  expect_error(ID(x = test_data, th.cons  = 0, contrast = "mean"))
  expect_error(ID(x = test_data, th.cons_lin = -10, contrast = "slope"))

  # The value for lambda is negative.
  expect_error(ID(x = test_data, lambda = -2))
  expect_error(ID(x = test_data, lambda.ic = -2))

  # The value for lambda is a positive real number instead of a positive integer.
  expect_warning(ID(x = test_data, lambda = 3.4))
  expect_warning(ID(x = test_data, lambda.ic = 5.4))
})
