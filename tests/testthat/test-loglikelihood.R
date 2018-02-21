library(vlad)
context("loglikelihood")
coeff <- c("(Intercept)"=-3.68, "Parsonnet"=0.077)

test_that("patients with different Parsonnet scores, RA=2, Steiner.etal (2000) p. 446", {
  expect_equal(round(loglikelihood(df=data.frame(as.integer(0), 0), coeff=coeff, RA=2), 3), -0.024)
  expect_equal(round(loglikelihood(df=data.frame(as.integer(0), 1), coeff=coeff, RA=2), 3), 0.67, tolerance=0.01)
  expect_equal(round(loglikelihood(df=data.frame(as.integer(50), 0), coeff=coeff, RA=2), 3), -0.43, tolerance=0.03)
  expect_equal(round(loglikelihood(df=data.frame(as.integer(50), 1), coeff=coeff, RA=2), 3), 0.26)
})

test_that("patients with different Parsonnet scores, RA=2, Steiner (2014) p. 234", {
  expect_equal(round(loglikelihood(df=data.frame(as.integer(0), 1), coeff=coeff, RA=2), 3), 0.669)
  expect_equal(round(loglikelihood(df=data.frame(as.integer(50), 0), coeff=coeff, RA=2), 3), -0.433)
})

test_that("patients with different Parsonnet scores, RA=1/2, Steiner (2014) p. 234", {
  expect_equal(round(loglikelihood(df=data.frame(as.integer(0), 0), coeff=coeff, RA=1/2), 3), 0.012)
  expect_equal(round(loglikelihood(df=data.frame(as.integer(0), 1), coeff=coeff, RA=1/2), 3), -0.681)
  expect_equal(round(loglikelihood(df=data.frame(as.integer(50), 0), coeff=coeff, RA=1/2), 3), 0.316)
  expect_equal(round(loglikelihood(df=data.frame(as.integer(50), 1), coeff=coeff, RA=1/2), 3), -0.377)
})
