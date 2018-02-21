library(vlad)
context("calceo")
coeff <- c("(Intercept)"=-3.68, "Parsonnet"=0.077)

test_that("patients with different Parsonnet scores, Steiner (2014) p.234", {
  expect_equal(round(calceo(df=data.frame(as.integer(0), 0), coeff=coeff)*-1, 3), -0.025)
  expect_equal(round(calceo(df=data.frame(as.integer(0), 1), coeff=coeff)*-1, 3), 0.975)
  expect_equal(round(calceo(df=data.frame(as.integer(50), 0), coeff=coeff)*-1, 3), -0.542)
  expect_equal(round(calceo(df=data.frame(as.integer(50), 1), coeff=coeff)*-1, 3), 0.458)
})


