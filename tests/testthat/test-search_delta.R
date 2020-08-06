context("search_delta")

data("cardiacsurgery", package = "spcadjust")

## preprocess data to 30 day mortality and subset data to
## phase I (In-control) and phase II (monitoring)
SALL <- cardiacsurgery %>% rename(s = Parsonnet) %>%
  mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
         phase = factor(ifelse(date < 2*365, "I", "II")))

tol <- 1e-6
## subset phase I (In-control)
SI <- filter(SALL, phase == "I") %>% select(s, y)

test_that("Search delta: ML", {
  expected_results <- .06304707
  works <- search_delta(SI$s, SI$y, type = "ML")
  expect_equal(works, expected_results, tolerance=tol)
})

test_that("Search delta: Pearson", {
  expected_results <- .07158091
  works <- search_delta(SI$s, SI$y, type = "Pearson")
  expect_equal(works, expected_results, tolerance=tol)
})
