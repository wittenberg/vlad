context("ell")

data("cardiacsurgery", package = "spcadjust")

## preprocess data to 30 day mortality and subset data to
## phase I (In-control) and phase II (monitoring)
SALL <- cardiacsurgery %>% mutate(s = Parsonnet) %>%
  mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
         phase = factor(ifelse(date < 2*365, "I", "II")))

## subset phase I (In-control)
SI <- filter(SALL, phase == "I") %>% select(s, y)

test_that("Maximum Likelihood", {
  expected_results <- -333.996
  dML <- search_delta(SI$s, SI$y, type = "ML")
  works <- ell(SI$s, SI$y, dML)
  expect_equal(as.numeric(works), expected_results, tolerance=1e-6)
})
