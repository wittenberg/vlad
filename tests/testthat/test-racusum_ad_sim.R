context("racusum_ad_sim")

library(vlad)
library(dplyr)
data("cardiacsurgery", package = "spcadjust")

set.seed(1234)
SALLI <- cardiacsurgery %>% mutate(s = Parsonnet) %>%
  mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
         phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
  filter(phase == "I") %>% select(s, y)

## estimate risk model, get relative frequences and probabilities
mod1 <- glm(y ~ s, data = SALLI, family = "binomial")
y <- SALLI$y
pi1 <- fitted.values(mod1)

## set up patient mix (risk model)
m <- 1e3
h <- 2
pmix <- data.frame(y, pi1, pi1)

test_that("Iterative search procedure I", {
  skip_on_cran()
  skip_if(SKIP==TRUE, "skip this test now")

  tol <- 0.3
  ## RA=2
  expected_results <- 424
  RLS <- sapply(1:m, racusum_ad_sim, pmix, h = h, RA = 2, RQ = 1, m = 50, type = "cond")
  works <- mean(RLS)
  expect_equal(works, expected_results, tolerance=tol)

  expected_results <- 443
  RLS <- sapply(1:m, racusum_ad_sim, pmix, h = 2, RA = 2, RQ = 1, m = 50, type = "cycl")
  works <- mean(RLS)
  expect_equal(works, expected_results, tolerance=tol)

  ## RA=1/2
  expected_results <- 576
  RLS <- sapply(1:m, racusum_ad_sim, pmix, h = 2, RA = 1/2, RQ = 1, m = 50, type = "cond")
  works <- mean(RLS)
  expect_equal(works, expected_results, tolerance=tol)

  expected_results <- 571
  RLS <- sapply(1:m, racusum_ad_sim, pmix, h = 2, RA = 1/2, RQ = 1, m = 50, type = "cycl")
  works <- mean(RLS)
  expect_equal(works, expected_results, tolerance=tol)
})
