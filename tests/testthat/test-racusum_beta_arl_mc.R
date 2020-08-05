context("racusum_arl_beta_arl_mc")

test_that("Different Markov Chain algorithms, detecting deterioration", {
  skip_on_cran()
  skip_if(SKIP == TRUE, "skip this test now")
  h <- 4.5
  RA <- 2
  g0 <- -3.6798
  g1 <- 0.0768*71
  shape1 <- 1
  shape2 <- 3
  r <- 600

  expected_results <- 4474.738
  MCtest <- list(
    racusum_beta_arl_mc(h=h, RA=RA, g0=g0, g1=g1, shape1=shape1, shape2=shape2, r=r, method=1),
    racusum_beta_arl_mc(h=h, RA=RA, g0=g0, g1=g1, shape1=shape1, shape2=shape2, r=r, method=2)
  )
  lapply(MCtest, function(x) expect_equal(x, expected_results, tolerance = 10^-6) )
})

test_that("Different Markov Chain algorithms, detecting improvement", {
  skip_on_cran()
  skip_if(SKIP == TRUE, "skip this test now")
  h <- 4.5
  RA <- 1/2
  g0 <- -3.6798
  g1 <- 0.0768*71
  shape1 <- 1
  shape2 <- 3
  r <- 600

  expected_results <- 5658.854
  MCtest <- list(
    racusum_beta_arl_mc(h=h, RA=RA, g0=g0, g1=g1, shape1=shape1, shape2=shape2, r=r, method=1),
    racusum_beta_arl_mc(h=h, RA=RA, g0=g0, g1=g1, shape1=shape1, shape2=shape2, r=r, method=2)
  )
  lapply(MCtest, function(x) expect_equal(x, expected_results, tolerance = 10^-6) )
})
