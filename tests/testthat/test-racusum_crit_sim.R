context("racusum_crit_sim")
R0 <- 1; RA <- 2
library("spcadjust")
data("cardiacsurgery")
cardiacsurgery <- dplyr::mutate(cardiacsurgery, phase=factor(ifelse(date < 2*365, "I", "II")))
S2 <- subset(cardiacsurgery, c(surgeon==2), c("phase", "Parsonnet", "status"))
S2I <- subset(S2, c(phase=="I"), c("Parsonnet", "status"))

df1 <- data.frame(Parsonnet=c(0L, 0L, 50L, 50L), status = c(0, 1, 0, 1))
coeff1 <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
L0 <- 1

test_that("Input parameter of function", {
  expect_error(racusum_crit_sim(L0 = 0, df1, coeff1))
})
#
test_that("Different input values for df", {
  dftest1 <- list(as.matrix(df1), NULL)
  lapply(dftest1, function(x) {
    expect_error(do.call(x, racusum_crit_sim(L0, df = x, coeff1)))})

  dftest2 <- list(data.frame(0L, 1, 1), data.frame(0L), data.frame(NA))
  lapply(dftest2, function(x) {
    expect_error(do.call(x, racusum_crit_sim(L0, df = x, coeff1)))})
})

test_that("Different input values for coeff", {
   coefftest <- list(coeff1[1], rep(1, 3), NULL, NA)
   lapply(coefftest, function(x) {
     expect_error(do.call(x, racusum_crit_sim(L0, df1, coeff = x)))})
})

test_that("Different input values for R0", {
  R0test <- list(-1, 0, "0", NA)
  lapply(R0test, function(x) {
    expect_error(do.call(x, racusum_crit_sim(L0, df1, coeff1, R0 = x)))})
})

test_that("Different input values for RA", {
  RAtest <- list(-1, 0, "0", NA)
  lapply(RAtest, function(x) {
    expect_error(do.call(x, racusum_crit_sim(L0, df1, coeff1, RA = x)))})
})

