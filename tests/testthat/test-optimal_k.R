context("optimal_k")

coeff <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
df <- data.frame(Parsonnet=c(0L, 0L, 50L, 50L),
                 status = c(0, 1, 0, 1))
QA <- 2
parsonnetscores <- df[, 1]

test_that("Output of optimal_k calculation", {
  expected_results <- 0.1255
  works <- round(optimal_k(QA, df, coeff, FALSE), 4)
  expect_equal(works, expected_results)
})

test_that("Different input values for coeff", {
  coefftest <- list(coeff[1], rep(1, 3), NULL, NA)
  lapply(coefftest, function(x) {
    expect_error(do.call(x, optimal_k(QA, df, coeff = x, FALSE)),
                 "Model coefficients 'coeff' must be a numeric vector with two elements")})
})

test_that("Different input values for QA", {
  QAtest <- list(-1, NA, 0)
  lapply(QAtest, function(x) {
    expect_error(do.call(x, optimal_k(QA = QAtest, df, coeff, FALSE)),
                 "QA must a positive numeric value")})
})
