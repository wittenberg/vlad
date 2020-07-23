context("eocusum_ad_sim")

df1 <- data.frame(Parsonnet=c(0L, 0L, 50L, 50L), status = c(0, 1, 0, 1))
coeff1 <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
coeff2 <- coeff1
k <- 0.01
r <- 1
h <- 1

# test_that("Input parameter of function", {
#   expect_error(eocusum_ad_sim(r = 0, k, h, df1, coeff1, coeff2))
#   expect_error(eocusum_ad_sim(r, k = -1, h, df1, coeff1, coeff2))
#   expect_error(eocusum_ad_sim(r, k, h = 0, df1, coeff1, coeff2))
# })
#
# test_that("Different input values for df", {
#   dftest1 <- list(as.matrix(df1), NULL)
#   lapply(dftest1, function(x) {
#     expect_error(do.call(x, eocusum_ad_sim(r, k, h, df = x, coeff1, coeff2)))})
#   })
#
#   dftest2 <- list(data.frame(0L, 1, 1), data.frame(0L), data.frame(NA))
#   lapply(dftest2, function(x) {
#     expect_error(do.call(x, eocusum_ad_sim(r, k, h, df = x, coeff1, coeff2)))
# })
#
#
# test_that("Different input values for coeff", {
#   coefftest <- list(coeff1[1], rep(1, 3), NULL, NA)
#   lapply(coefftest, function(x) {
#     expect_error(do.call(x, eocusum_ad_sim(r, k, h, df1, coeff = x, coeff2)))})
# })
#
# test_that("Different input values for coeff2", {
#   coefftest <- list(coeff2[1], rep(1, 3), NULL, NA)
#   lapply(coefftest, function(x) {
#     expect_error(do.call(x, eocusum_ad_sim(r, k, h, df1, coeff1, coeff2 = x)))})
# })
#
# test_that("Input parameter QS", {
#   QStest <- list(-1, 0)
#   lapply(QStest, function(x) {
#     expect_error(do.call(x, eocusum_ad_sim(r, k, h, df1, coeff1, coeff2, QS = x)))})
#   expect_error(eocusum_ad_sim(r, k, h, df1, coeff1, coeff2, QS = 1/2, side = "low"),
#                fixed = TRUE)
#   expect_error(eocusum_ad_sim(r, k, h, df1, coeff1, coeff2, QS = 2, side = "up"),
#                fixed = TRUE)
# })
#
# test_that("Input value for side", {
#   sidetest <- "A"
#   expect_error(eocusum_ad_sim(r, k, h, df1, coeff1, coeff2, side = sidetest), fixed = TRUE)
# })

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

  # type = cond
  kopt_det <- optimal_k(pmix, RA=2)
  RLS <- sapply(1:m, eocusum_ad_sim, pmix, k=kopt_det, h=2, RQ = RQ, type = "cond", side = "low")
  expect_equal(mean(RLS), 200, tolerance=tol)

  kopt_imp <- optimal_k(pmix, RA=1/2)
  RLS <- sapply(1:m, eocusum_ad_sim, pmix, k=kopt_imp, h=2, RQ = RQ, type = "cond", side = "up")
  expect_equal(mean(RLS), 265, tolerance=tol)

  # type = cycl
  # kopt_det <- optimal_k(pmix, RA=2)
  # RLS <- sapply(1:m, eocusum_ad_sim, pmix, k=kopt_det, h=2, RQ = RQ, type = "cycl", side = "low")
  # expect_equal(mean(RLS), 197, tolerance=tol)

  kopt_imp <- optimal_k(pmix, RA=1/2)
  RLS <- sapply(1:m, eocusum_ad_sim, pmix, k=kopt_imp, h=2, RQ = RQ, type = "cycl", side = "up")
  expect_equal(mean(RLS), 160, tolerance=tol)
})
