context("eocusum_arl_sim")

df1 <- data.frame(Parsonnet=c(0L, 0L, 50L, 50L), status = c(0, 1, 0, 1))
coeff1 <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
k <- 0.01
r <- 1
h <- 1

test_that("Input parameter of function", {
  expect_error(eocusum_arl_sim(r = 0, k, h, df1, coeff1))
  expect_error(eocusum_arl_sim(r, k = -1, h, df1, coeff1))
  expect_error(eocusum_arl_sim(r, k, h = 0, df1, coeff1))
})

test_that("Different input values for df", {
  dftest1 <- list(as.matrix(df1), NULL)
  lapply(dftest1, function(x) {
    expect_error(do.call(x, eocusum_arl_sim(r, k, h, df = x, coeff1)))})

  dftest2 <- list(data.frame(0L, 1, 1), data.frame(0L), data.frame(NA))
  lapply(dftest2, function(x) {
    expect_error(do.call(x, eocusum_arl_sim(r, k, h, df = x, coeff1)))})

  dftest3 <- list(data.frame(0, 1), data.frame("0", 1), data.frame(NA, 1))
  lapply(dftest3, function(x) {
    expect_error(do.call(x, eocusum_arl_sim(r, k, h, df = x, coeff1)))})

  dftest4 <- list(data.frame(0L, 1L), data.frame(0L, "1L"), data.frame(0L, NA))
  lapply(dftest4, function(x) {
    expect_error(do.call(x, eocusum_arl_sim(r, k, h, df = x, coeff1)))})
})

test_that("Different input values for coeff", {
  coefftest <- list(coeff1[1], rep(1, 3), NULL, NA)
  lapply(coefftest, function(x) {
    expect_error(do.call(x, eocusum_arl_sim(r, k, h, df1, coeff = x)))})
})

test_that("Different input values for yemp", {
  expect_error(eocusum_arl_sim(r, k, h, df1, coeff1, yemp = as.character(TRUE)))
  expect_error(eocusum_arl_sim(r, k, h, df1, coeff1, yemp = as.numeric(TRUE)))
  expect_error(eocusum_arl_sim(r, k, h, df1, coeff1, yemp = NA))
})

test_that("Input value for side", {
  sidetest <- "A"
  expect_error(eocusum_arl_sim(r, k, h, df1, coeff1, side = sidetest), fixed=TRUE)
})

test_that("Iterative search procedure I", {
  skip_on_cran()
  skip_if(SKIP==TRUE, "skip this test now")

  data("cardiacsurgery", package = "spcadjust")
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
  h <- 2
  m <- 1e3

  set.seed(1234)
  RQ <- 1
  expected_results <- 1000
  m <- 1e4
  tol <- 0.3

  # yemp = FALSE
  kopt_det <- optimal_k(pmix, RA=2)
  RLS <- sapply(1:m, eocusum_arl_sim, pmix, k=kopt_det, h=2, RQ = RQ, yemp = FALSE, side = "low")
  expect_equal(mean(RLS), 210, tolerance=tol)

  kopt_imp <- optimal_k(pmix, RA=1/2)
  RLS <- sapply(1:m, eocusum_arl_sim, pmix, k=kopt_imp, h=2, RQ = RQ, yemp = FALSE, side = "up")
  expect_equal(mean(RLS), 304, tolerance=tol)

  # yemp = TRUE
  kopt_det <- optimal_k(pmix, RA=2)
  RLS <- sapply(1:m, eocusum_arl_sim, pmix, k=kopt_det, h=2, RQ = RQ, yemp = TRUE, side = "low")
  expect_equal(mean(RLS), 197, tolerance=tol)

  kopt_imp <- optimal_k(pmix, RA=1/2)
  RLS <- sapply(1:m, eocusum_arl_sim, pmix, k=kopt_imp, h=2, RQ = RQ, yemp = TRUE, side = "up")
  expect_equal(mean(RLS), 286, tolerance=tol)
})
















