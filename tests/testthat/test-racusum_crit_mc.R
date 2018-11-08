context("racusum_arl_h_mc")
library("dplyr")
data("cardiacsurgery", package = "spcadjust")

## preprocess data to 30 day mortality and subset phase I (In-control) of surgeons 2
SALLI <- cardiacsurgery %>% rename(s = Parsonnet) %>%
  mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
         phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
  filter(phase == "I") %>% select(s, y)

## estimate risk model, get relative frequences and probabilities
mod1 <- glm(y ~ s, data = SALLI, family = "binomial")
fi  <- as.numeric(table(SALLI$s) / length(SALLI$s))
usi <- sort(unique(SALLI$s))
pi <- predict(mod1, newdata = data.frame(s = usi), type = "response")

## set up patient mix
pmix  <- data.frame(fi, pi, pi)

RQ <- 1

## control limit for detecting deterioration R1=2:

test_that("Iterative search procedure I, deterioration R1 = 2", {
  skip_on_cran()
  skip_if(SKIP == TRUE, "skip this test now")
  RA <- 2
  tol <- 0.01
  L0 <- 4*10^4
  ## paired rounding
  expected_results <- racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "p", method = "Toep", verbose = TRUE)
  MCtest <- list(
    racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "p", method = "ToepInv", verbose = TRUE),
    racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "p", method = "BE", verbose = TRUE)
  )
  lapply(MCtest, function(x) expect_equal(x, expected_results, tolerance = tol) )
  ## simple rounding
  expected_results <- racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "s", method = "Toep", verbose = TRUE)
  MCtest <- list(
    racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "s", method = "ToepInv", verbose = TRUE),
    racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "s", method = "BE", verbose = TRUE)
  )
  lapply(MCtest, function(x) expect_equal(x, expected_results, tolerance = tol) )
  ## paired rounding, extra cases
  L0 <- 10^5
  expected_results <- racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "p", method = "Toep", verbose = TRUE)
  MCtest <- list(
    racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "p", method = "ToepInv", verbose = TRUE),
    racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "p", method = "BE", verbose = TRUE)
  )
  lapply(MCtest, function(x) expect_equal(x, expected_results, tolerance = tol) )
  ## simple rounding, extra cases
  expected_results <- racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "s", method = "Toep", verbose = TRUE)
  MCtest <- list(
    racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "s", method = "ToepInv", verbose = TRUE),
    racusum_crit_mc(pmix = pmix, L0 = L0, RA = RA, RQ = RQ, rounding = "s", method = "BE", verbose = TRUE)
  )
  lapply(MCtest, function(x) expect_equal(x, expected_results, tolerance = tol) )
})
