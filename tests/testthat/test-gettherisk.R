library(vlad)
context("gettherisk")
coeff <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)

test_that("Different Parsonnet scores, Steiner et al. (2000) p. 445", {
  expect_equal(round(gettherisk(0, coeff), 3), 0.025)
  expect_equal(round(gettherisk(71, coeff), 3), 0.86, tolerance = 0.03)
})

test_that("Different Parsonnet scores, Steiner (2014) p. 234", {
  expect_equal(round(gettherisk(0, coeff), 3), 0.025)
  expect_equal(round(gettherisk(50, coeff), 3), 0.542)
})
