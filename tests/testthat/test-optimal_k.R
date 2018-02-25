library(vlad)
context("optimal_k")
coeff <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
df <- data.frame(Parsonnet=c(0L, 0L, 50L, 50L),
                 status = c(0, 1, 0, 1))
QA <- 2

test_that("Output of optimal_k calculation", {
  expected_results <- 0.1255
  works <- round(optimal_k(QA, parsonnetscores = df[, 1], coeff), 4)
  expect_equal(works, expected_results)
})

test_that("Different input values for coeff", {
  coeff3 <- list(coeff[1], rep(1, 3), NULL)
  lapply(coeff3, function(x) {
    expect_error(do.call(x, optimal_k(QA, df[, 1], coeff = coeff3)),
                 "model coefficients \"coeff\" must a numeric vector with two elements")})
})
