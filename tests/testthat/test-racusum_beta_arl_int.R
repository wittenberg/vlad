context("racusum_beta_arl_int")

test_that("Different integration algorithms, detecting deterioration", {
  skip_on_cran()
  skip_if(SKIP == TRUE, "skip this test now")
  h <- 4.5
  N <- 70
  RQ <- 1
  g0 <- -3.6798
  g1 <- 0.0768*71
  shape1 <- 1
  shape2 <- 3
  tol <- 10^-6
  expect_equal(racusum_beta_arl_int(h=h, N=N, RA=2, RQ=RQ, g0=g0, g1=g1, shape1=shape1, shape2=shape2, pw=TRUE), 4485.203, tolerance=tol)
  expect_equal(racusum_beta_arl_int(h=h, N=N, RA=2, RQ=RQ, g0=g0, g1=g1, shape1=shape1, shape2=shape2, pw=FALSE), 4561.862, tolerance=tol)
  expect_equal(racusum_beta_arl_int(h=h, N=N, RA=1/2, RQ=RQ, g0=g0, g1=g1, shape1=shape1, shape2=shape2, pw=TRUE), 5731.772, tolerance=tol)
  expect_equal(racusum_beta_arl_int(h=h, N=N, RA=1/2, RQ=RQ, g0=g0, g1=g1, shape1=shape1, shape2=shape2, pw=FALSE), 5728.431, tolerance=tol)
})
