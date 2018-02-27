context("llr_score")

coeff <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
df <- data.frame(Parsonnet = c(0L, 0L, 50L, 50L),
                status = c(0, 1, 0, 1))

RA <- 2
test_that("patients with different Parsonnet scores, RA = 2, Steiner et al. (2000) p. 446", {
  expected_results <- list(-0.024, 0.67, -0.43, 0.26)
  works <- lapply(seq_along(df$Parsonnet), function(i) round(llr_score(df = df[i, ], coeff = coeff, RA = RA), 3))
  expect_equal(works, expected_results, tolerance = 0.03)
})

test_that("patients with different Parsonnet scores, RA = 2, Steiner (2014) p. 234", {
  expected_results <- list(-0.024, 0.669, -0.433, 0.26)
  works <- lapply(seq_along(df$Parsonnet), function(i) round(llr_score(df = df[i, ], coeff = coeff, RA = RA), 3))
  expect_equal(works, expected_results)
})

RA <- 1/2
test_that("patients with different Parsonnet scores, RA = 1/2, Steiner (2014) p. 234", {
  expected_results <- list(0.012, -0.681, 0.316, -0.377)
  works <- lapply(seq_along(df$Parsonnet), function(i) round(llr_score(df = df[i, ], coeff = coeff, RA = RA), 3))
  expect_equal(works, expected_results)
})

test_that("Different input values for df", {
  dftest1 <- list(as.matrix(df), NULL)
  lapply(dftest1, function(x) {
    expect_error(do.call(x, llr_score(df = x, coeff)),
                 "Provide a dataframe for argument 'df'")})

  dftest2 <- list(data.frame(0L, 1, 1), data.frame(0L), data.frame(NA))
  lapply(dftest2, function(x) {
    expect_error(do.call(x, llr_score(df = x, coeff)),
                 "Provide a dataframe with two columns for argument 'df'")})

  dftest3 <- list(data.frame(0, 1), data.frame("0", 1), data.frame(NA, 1))
  lapply(dftest3, function(x) {
    expect_error(do.call(x, llr_score(df = x, coeff)),
                 "First column of dataframe must be of type integer")})

  dftest4 <- list(data.frame(0L, 1L), data.frame(0L, "1L"), data.frame(0L, NA))
  lapply(dftest4, function(x) {
    expect_error(do.call(x, llr_score(df = x, coeff)),
                 "Second column of dataframe must be of type numeric")})
})

test_that("Different input values for coeff", {
  coeff3 <- list(coeff[1], rep(1, 3), NULL, NA)
  lapply(coeff3, function(x) {
    expect_error(do.call(x, llr_score(df, coeff = x)),
                 "Model coefficients 'coeff' must be a numeric vector with two elements")})
})

test_that("Different input values for R0", {
  R0test <- list(-1, 0, "0", NA)
  lapply(R0test, function(x) {
    expect_error(do.call(x, llr_score(df, coeff, R0 = x)),
                 "Odds ratio of death under the null hypotheses 'R0' must a positive numeric value")})
})

test_that("Different input values for RA", {
  R0test <- list(-1, 0, "0", NA)
  lapply(R0test, function(x) {
    expect_error(do.call(x, llr_score(df, coeff, RA = x)),
                 "Odds ratio of death under the alternative hypotheses 'RA' must a positive numeric value")})
})

test_that("Different input values for yemp", {
  expect_warning(llr_score(df, coeff, yemp = as.character(TRUE)),
                 "Argument 'yemp' must be logical using TRUE as default value")
  expect_warning(llr_score(df, coeff, yemp = as.numeric(TRUE)),
                 "Argument 'yemp' must be logical using TRUE as default value")
  expect_warning(llr_score(df, coeff, yemp = NA),
                 "Argument 'yemp' must be logical using TRUE as default value")
})
