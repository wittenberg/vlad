library(vlad)
context("eocusum_arloc_sim")

df1 <- data.frame(Parsonnet=c(0L, 0L, 50L, 50L), status = c(0, 1, 0, 1))
coeff1 <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
coeff2 <- coeff1
k <- 0.01
r <- 1
h <- 1

test_that("Input parameter of function", {
  expect_error(eocusum_arloc_sim(r = 0, k, h, df1, coeff1, coeff2),
               "number of simulation runs r must a positive integer")
  expect_error(eocusum_arloc_sim(r, k = -1, h, df1, coeff1, coeff2),
               "reference value k must a positive numeric value")
  expect_error(eocusum_arloc_sim(r, k, h = 0, df1, coeff1, coeff2),
               "control limit h must a positive numeric value")
})

test_that("Different input values for df", {
  expect_error(eocusum_arloc_sim(r, k, h, df = NULL, coeff1),
               "provide a dataframe with two columns for argument \"df\"")
  expect_error(eocusum_arloc_sim(r, k, h, df = data.frame(0L, as.character(1)), coeff1, coeff2),
               "second column of dataframe must be of type numeric")
  expect_error(eocusum_arloc_sim(r, k, h, df = data.frame(0L, as.integer(1)), coeff1, coeff2),
               "second column of dataframe must be of type numeric")
  expect_error(eocusum_arloc_sim(r, k, h, df = data.frame(as.character(0L), 1), coeff1, coeff2),
               "first column of dataframe must be of type integer")
})

test_that("Different input values for coeff1", {
  expect_error(eocusum_arloc_sim(r, k, h, df1, coeff = coeff1[1], coeff2),
               "model coefficients \"coeff\" must a numeric vector with two elements")
  expect_error(eocusum_arloc_sim(r, k, h, df1, coeff = rep(1, 3), coeff2),
               "model coefficients \"coeff\" must a numeric vector with two elements")
  expect_error(eocusum_arloc_sim(r, k, h, df1, coeff = NULL, coeff2),
               "model coefficients \"coeff\" must a numeric vector with two elements")
})

test_that("Different input values for coeff2", {
  expect_error(eocusum_arloc_sim(r, k, h, df1, coeff1, coeff2 = coeff2[1]),
               "model coefficients \"coeff2\" must a numeric vector with two elements")
  expect_error(eocusum_arloc_sim(r, k, h, df1, coeff1, coeff2 = rep(1, 3)),
               "model coefficients \"coeff2\" must a numeric vector with two elements")
  expect_error(eocusum_arloc_sim(r, k, h, df1, coeff1, coeff2 = NULL),
               "model coefficients \"coeff2\" must a numeric vector with two elements")
})

# test_that("Different input values for yemp", {
#   # expect_warning(eocusum_arloc_sim(r, k, h, df1, coeff1, coeff2, yemp = as.character(TRUE)),
#   #                "argument \"yemp\" must be logical using TRUE as default value")
#   expect_warning(eocusum_arloc_sim(r, k, h, df1, coeff1, coeff2, yemp = as.numeric(TRUE)),
#                  "argument \"yemp\" must be logical using TRUE as default value")
# })

# test_that("Different input values for side", {
#   expect_warning(eocusum_arl_sim(r, kopt, h, df1, coeff1, side=2),
#                  "no valid input, using side=low (deterioration) as default")
# })
