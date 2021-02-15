context("racusum_crit_mc")

pmix <- data.frame(cbind(c(0.1, 0.5, 0.4), c(0.1, 0.12, 0.14), c(0.1, 0.12, 0.14)))
RA <- 2
RQ <- 1
h <- 2.9
scaling <- 600
rounding <- "p"
method <- "Toep"
L0 <- 100

test_that("Different input values for pmix", {
  pmix2 <- data.frame(cbind(c(0.2, 0.5, 0.4), c(0.1, 0.12, 0.14), c(0.1, 0.12, 0.14)))
  expect_error(racusum_crit_mc(pmix = pmix2, L0, RA, RQ, h))
})

# test_that("Different input values for RA", {
#   RAtest <- list(-1, 0, "0", NA)
#   lapply(RAtest, function(x) {
#     expect_error(do.call(x, racusum_crit_mc(pmix, L0, RQ, scaling, rounding, method, RA = x)))})
# })
#
test_that("Different input values for RQ", {
  R0test <- list(-1, 0, "0", NA)
  lapply(R0test, function(x) {
    expect_error(do.call(x, racusum_crit_mc(pmix, L0, RA, scaling, rounding, method, RQ = x)))})
})

# test_that("Input parameter for L0", {
#   expect_error(racusum_crit_mc(pmix, L0 = 0, RA, RQ, scaling, rounding, method))
# })

test_that("Different input values for scaling", {
  scatest <- list(-1, 0, "0")
  lapply(scatest, function(x) {
    expect_error(do.call(x, racusum_crit_mc(pmix, L0, RA, RQ, rounding, method, scaling = x)))})
})

test_that("Different Markov Chain algorithms, detecting deterioration", {
  skip_on_cran()
  skip_if(SKIP == TRUE, "skip this test now")
  data("cardiacsurgery", package = "spcadjust")
  ## preprocess data to 30 day mortality and subset phase I (In-control) of surgeons 2
  SALLI <- cardiacsurgery %>% mutate(s = Parsonnet) %>%
    mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
           phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
    filter(phase == "I") %>% select(s, y)
  # set up parameters
  mod1 <- glm(y ~ s, data = SALLI, family = "binomial")
  s <- sort(unique(SALLI$s))
  fi  <- as.numeric(table(SALLI$s)/length(SALLI$s))
  pi <- predict(mod1, newdata = data.frame(s), type = "response")
  pmix  <- data.frame(fi, pi, pi)

  L0 <- 370
  sca <- 600
  RA <- 2
  RQ <- 1

  works <- 1.8592
  ## paired rounding
  MCtest <- list(
    racusum_crit_mc(L0=L0, pmix = pmix, RA = RA, RQ = RQ, scaling = sca, rounding = "p", method = "Toep", verbose=FALSE),
    racusum_crit_mc(L0=L0, pmix = pmix, RA = RA, RQ = RQ, scaling = sca, rounding = "p", method = "Toep", verbose=TRUE)
  )
  lapply(MCtest, function(x) expect_equal(x, works, tolerance = 10^-4) )
})
