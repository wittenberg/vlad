context("QQ")

data("cardiacsurgery", package = "spcadjust")

## preprocess data to 30 day mortality and subset data to
## phase I (In-control) and phase II (monitoring)
SALL <- cardiacsurgery %>% mutate(s = Parsonnet) %>%
  mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
         phase = factor(ifelse(date < 2*365, "I", "II")))

## subset phase I (In-control)
SI <- filter(SALL, phase == "I") %>% select(s, y)

test_that("Pearson measure", {
  expected_results <- 61.01746
  dQQ <- search_delta(SI$s, SI$y, type = "Pearson")
  works <- QQ(SI$s, SI$y, dQQ)
  expect_equal(works, expected_results, tolerance=1e-6)
})
