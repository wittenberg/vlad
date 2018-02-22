library(vlad)
context("loglikelihood")
coeff <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
df <- data.frame(Parsonnet=c(0L, 0L, 50L, 50L),
                status = c(0, 1, 0, 1))

RA <- 2
test_that("patients with different Parsonnet scores, RA = 2, Steiner et al. (2000) p. 446", {
  expected_results <- list(-0.024, 0.67, -0.43, 0.26)
  works <- lapply(1:nrow(df), function(i) round(loglikelihood(df=df[i, ], coeff = coeff, RA = RA), 3))
  expect_equal(works, expected_results, tolerance = 0.03)
})

test_that("patients with different Parsonnet scores, RA = 2, Steiner (2014) p. 234", {
  expected_results <- list(-0.024, 0.669, -0.433, 0.26)
  works <- lapply(1:nrow(df), function(i) round(loglikelihood(df=df[i, ], coeff = coeff, RA = RA), 3))
  expect_equal(works, expected_results)
})

RA <- 1/2
test_that("patients with different Parsonnet scores, RA = 1/2, Steiner (2014) p. 234", {
  expected_results <- list(0.012, -0.681, 0.316, -0.377)
  works <- lapply(1:nrow(df), function(i) round(loglikelihood(df=df[i, ], coeff = coeff, RA = RA), 3))
  expect_equal(works, expected_results)
})






