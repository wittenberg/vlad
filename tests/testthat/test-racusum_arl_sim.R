context("racusum_arl_sim")

df1 <- data.frame(Parsonnet=c(0L, 0L, 50L, 50L), status = c(0, 1, 0, 1))
coeff1 <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
r <- 1
h <- 1

test_that("Input parameter of function", {
  expect_error(racusum_arl_sim(r = 0, coeff1, h, df1))
  expect_error(racusum_arl_sim(r, coeff1, h = 0, df1))
})

dftest2 <- list(data.frame(0L, 1, 1), data.frame(0L), data.frame(NA))
  lapply(dftest2, function(x) {
  expect_error(do.call(x, racusum_arl_sim(r, coeff1, h, df = x)))})

dftest3 <- list(data.frame(0, 1), data.frame("0", 1), data.frame(NA, 1))
  lapply(dftest3, function(x) {
  expect_error(do.call(x, racusum_arl_sim(r, coeff1, h, df = x)))})

dftest4 <- list(data.frame(0L, 1L), data.frame(0L, "1L"), data.frame(0L, NA))
  lapply(dftest4, function(x) {
  expect_error(do.call(x, racusum_arl_sim(r, coeff1, h, df = x)))})

test_that("Different input values for coeff", {
  coefftest <- list(coeff1[1], rep(1, 3), NULL, NA)
  lapply(coefftest, function(x) {
    expect_error(do.call(x, racusum_arl_sim(r, coeff = x, h, df1)))})
})

test_that("Different input values for RA", {
  RAtest <- list(-1, 0, "0", NA)
  lapply(RAtest, function(x) {
    expect_error(do.call(x, racusum_arl_sim(r, coeff1, h, df1, RA = x)))})
})

test_that("Different input values for yemp", {
  expect_error(racusum_arl_sim(r, coeff1, h, df1, yemp = as.character(TRUE)))
  expect_error(racusum_arl_sim(r, coeff1, h, df1, yemp = as.numeric(TRUE)))
  expect_error(racusum_arl_sim(r, coeff1, h, df1, yemp = NA))
})

test_that("Iterative search procedure I", {
  skip_on_cran()
  skip_if(SKIP==TRUE, "skip this test now")
  ## preprocess data to 30 day mortality and subset phase I (In-control) of surgeons 2
  SALLI <- cardiacsurgery %>% mutate(s = Parsonnet) %>%
    mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
        phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
    filter(phase == "I") %>% select(s, y)

  ## estimate risk model, get relative frequences and probabilities
  mod1 <- glm(y ~ s, data = SALLI, family = "binomial")
  y <- SALLI$y
  pi1 <- fitted.values(mod1)

  ## set up patient mix (risk model)
  pmix <- data.frame(y, pi1, pi1)
  h <- 2.75599

  set.seed(1234)
  ## RA=2
  expected_results <- 1000
  m <- 1e4
  RLS <- sapply(1:m, racusum_arl_sim, h=h, pmix=pmix, RA=2, yemp=TRUE)
  works <- mean(RLS)
  expect_equal(works, expected_results, tolerance=0.3)

  expected_results <- 1150
  m <- 1e4
  RLS <- sapply(1:m, racusum_arl_sim, h=h, pmix=pmix, RA=2, yemp=FALSE)
  works <- mean(RLS)
  expect_equal(works, expected_results, tolerance=0.3)

  ## RA=1/2
  expected_results <- 1450
  m <- 1e4
  RLS <- sapply(1:m, racusum_arl_sim, h=h, pmix=pmix, RA=1/2, yemp=TRUE)
  works <- mean(RLS)
  expect_equal(works, expected_results, tolerance=0.3)

  expected_results <- 1600
  m <- 1e4
  RLS <- sapply(1:m, racusum_arl_sim, h=h, pmix=pmix, RA=1/2, yemp=FALSE)
  works <- mean(RLS)
  expect_equal(works, expected_results, tolerance=0.3)
})
