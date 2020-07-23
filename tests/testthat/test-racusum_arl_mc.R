context("racusum_arl_mc")

pmix <- data.frame(cbind(c(0.1, 0.5, 0.4), c(0.1, 0.12, 0.14), c(0.1, 0.12, 0.14)))
RA <- 2
RQ <- 1
h <- 2.9
scaling <- 600
rounding <- "p"
method <- "Toep"

test_that("Different input values for RA", {
  RAtest <- list(-1, 0, "0", NA)
  lapply(RAtest, function(x) {
    expect_error(do.call(x, racusum_arl_mc(pmix, RQ, h, scaling, rounding, method, RA = x)))
    })
})

test_that("Different input values for RQ", {
  R0test <- list(-1, 0, "0", NA)
  lapply(R0test, function(x) {
    expect_error(do.call(x, racusum_arl_mc(pmix, RA, h, scaling, rounding, method, RQ = x)))
    })
})

test_that("Input parameter for h", {
  expect_error(racusum_arl_mc(pmix, RA, RQ, h = 0, scaling, rounding, method))
})

test_that("Different input values for scaling", {
  scatest <- list(-1, 0, "O")
  lapply(scatest, function(x) {
    expect_error(do.call(x, racusum_arl_mc(pmix, RA, RQ, h, rounding, method, scaling = x)))}
    )
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

  h <- 4.655
  sca <- 600

  works <- racusum_arl_mc(pmix = pmix, h = h, RA = 2, RQ = 1, scaling = sca, rounding = "s", method = "Toep")
  ## simple rounding
  MCtest <- list(
    racusum_arl_mc(pmix = pmix, h = h, RA = 2, RQ = 1, scaling = sca, rounding = "s", method = "ToepInv"),
    racusum_arl_mc(pmix = pmix, h = h, RA = 2, RQ = 1, scaling = sca, rounding = "s", method = "BE")
  )
  lapply(MCtest, function(x) expect_equal(x, works, tolerance = 10^-6) )

  ## paired rounding
  works <- racusum_arl_mc(pmix = pmix, h = h, RA = 2, RQ = 1, scaling = sca, rounding = "p", method = "Toep")
  MCtest <- list(
    racusum_arl_mc(pmix = pmix, h = h, RA = 2, RQ = 1, scaling = sca, rounding = "p", method = "ToepInv"),
    racusum_arl_mc(pmix = pmix, h = h, RA = 2, RQ = 1, scaling = sca, rounding = "p", method = "BE")
  )
  lapply(MCtest, function(x) expect_equal(x, works, tolerance = 10^-6) )
})

test_that("Different Markov Chain algorithms, detecting improvement", {
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
  h <- 4.655
  sca <- 600

  works <- racusum_arl_mc(pmix = pmix, h = h, RA = 1/2, RQ = 1, scaling = sca, rounding = "s", method = "Toep")
  ## simple rounding
  MCtest <- list(
    racusum_arl_mc(pmix = pmix, h = h, RA = 1/2, RQ = 1, scaling = sca, rounding = "s", method = "ToepInv"),
    racusum_arl_mc(pmix = pmix, h = h, RA = 1/2, RQ = 1, scaling = sca, rounding = "s", method = "BE")
  )
  lapply(MCtest, function(x) expect_equal(x, works, tolerance = 10^-6) )

  works <- racusum_arl_mc(pmix = pmix, h = h, RA = 1/2, RQ = 1, scaling = sca, rounding = "p", method = "Toep")
  ## paired rounding
  MCtest <- list(
    racusum_arl_mc(pmix = pmix, h = h, RA = 1/2, RQ = 1, scaling = sca, rounding = "p", method = "ToepInv"),
    racusum_arl_mc(pmix = pmix, h = h, RA = 1/2, RQ = 1, scaling = sca, rounding = "p", method = "BE")
  )
  lapply(MCtest, function(x) expect_equal(x, works, tolerance = 10^-6) )
})
