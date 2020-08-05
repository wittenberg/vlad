context("racusum_arl_beta_crit_mc")

test_that("Different Markov Chain algorithms, detecting deterioration", {
  skip_on_cran()
  skip_if(SKIP == TRUE, "skip this test now")
  L0 <- 500
  g0 <- -3.6798
  g1 <- 0.0768*71
  shape1 <- 1
  shape2 <- 3
  r <- 600
  tol <- 10^-6

  expected_results <- 2.5178
  MCtest <- list(
    racusum_beta_crit_mc(L0=L0, RA=2, g0=g0, g1=g1, shape1=shape1, shape2=shape2, r=r, method=1, verbose=TRUE),
    racusum_beta_crit_mc(L0=L0, RA=2, g0=g0, g1=g1, shape1=shape1, shape2=shape2, r=r, method=1, verbose=FALSE)
  )
  lapply(MCtest, function(x) expect_equal(x, expected_results, tolerance = tol) )
})
