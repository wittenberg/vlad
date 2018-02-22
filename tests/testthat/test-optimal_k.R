library(vlad)
context("optimal_k")
coeff <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
df <- data.frame(Parsonnet=c(0L, 0L, 50L, 50L),
                 status = c(0, 1, 0, 1))
RA <- 2

test_that("Output of optimal_k calculation", {
  expected_results <- 0.1255
  works <- round(optimal_k(parsonnetscores=df[, 1], coeff = coeff, QA = RA), 4)
  expect_equal(works, expected_results)
})
