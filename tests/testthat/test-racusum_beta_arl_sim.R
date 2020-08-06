context("racusum_beta_arl_sim")

set.seed(1234)
h <- 4.5
RQ <- 1
maxS <- 71
g0 <- -3.6798
g1 <- 0.0768
shape1 <- 1
shape2 <- 3
tol <- 10^-6

test_that("Different input values for RA", {
  RAtest <- list(-1, 0, "0", NA)
  lapply(RAtest, function(x) {
    expect_error(do.call(x, racusum_beta_arl_sim, h=h, RQ=RQ, coeff=c(g0,g1), shape1=shape1, shape2=shape2, rs=maxS, RA=x))})
})

test_that("Different simulation algorithms, detecting deterioration", {
  skip_on_cran()
  skip_if(SKIP == TRUE, "skip this test now")

  ## RA=2
  expected_results <- 4568
  m <- 1e4
  RLS <- sapply(1:m, racusum_beta_arl_sim, h=h, RA=2, RQ=RQ, coeff=c(g0,g1), shape1=shape1, shape2=shape2, rs=maxS)
  works <- mean(RLS)
  expect_equal(works, expected_results, tolerance=0.3)

  ## RA=1/2
  expected_results <- 5644
  m <- 1e4
  RLS <- sapply(1:m, racusum_beta_arl_sim, h=h, RA=1/2, RQ=RQ, coeff=c(g0,g1), shape1=shape1, shape2=shape2, rs=maxS)
  works <- mean(RLS)
  expect_equal(works, expected_results, tolerance=0.3)
})
