library(vlad)
context("calceo")
coeff <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
df <- data.frame(Parsonnet=c(0L, 0L, 50L, 50L),
                 status = c(0, 1, 0, 1))

test_that("Output results for different Parsonnet scores, Steiner (2014) p.234", {
  expected_results <- list(-0.025, 0.975, -0.542, 0.458)
  works <- lapply(seq_along(df$Parsonnet), function(i) round(calceo(df=df[i, ], coeff = coeff)*-1, 3))
  expect_equal(works, expected_results)
})

test_that("Different input values for coeff", {
  coeff3 <- list(coeff[1], rep(1, 3), NULL)
  lapply(coeff3, function(x) {
    expect_error(do.call(x, calceo(df, coeff = coeff3)),
                 "model coefficients 'coeff' must a numeric vector with two elements")})
})

test_that("Different input values for df", {
  expect_error(calceo(df = NULL, NULL),
               "provide a dataframe with two columns for argument 'df'")
  expect_error(calceo(df = NULL, coeff),
               "provide a dataframe with two columns for argument 'df'")
  expect_error(calceo(df = data.frame(0L, as.character(1)), coeff),
               "second column of dataframe must be of type numeric")
  expect_error(calceo(df = data.frame(0L, as.integer(1)), coeff),
               "second column of dataframe must be of type numeric")
  expect_error(calceo(df = data.frame(as.character(0L), 1), coeff),
               "first column of dataframe must be of type integer")
})

test_that("Different input values for yemp", {
  expect_warning(calceo(df, coeff, yemp = as.character(TRUE)),
                 "argument 'yemp' must be logical using TRUE as default value")
  expect_warning(calceo(df, coeff, yemp = as.numeric(TRUE)),
                 "argument 'yemp' must be logical using TRUE as default value")
})
