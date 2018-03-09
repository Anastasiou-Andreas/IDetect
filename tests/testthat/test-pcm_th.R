context("The change-points for piecewise-constant mean through thresholding")

# In this file we test the pcm_th function

test_that("The correct change-points are given", {
  set.seed(100)
  no.cpt <- rnorm(100)
  one.cpt <- rnorm(100) + c(rep(-10, 50), rep(20, 50))
  two.cpts <- rnorm(500) + c(rep(0, 200), rep(10, 150), rep(-10, 150))
  expect_equal(pcm_th(no.cpt), numeric(0))
  expect_equal(pcm_th(one.cpt), 50)
  expect_equal(pcm_th(two.cpts), c(200, 350))
})
