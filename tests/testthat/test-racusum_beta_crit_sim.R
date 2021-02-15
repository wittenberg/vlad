context("racusum_beta_crit_sim")

L0 <- 500
RQ <- 1
maxS <- 71
g0 <- -3.6798
g1 <- 0.0768
shape1 <- 1
shape2 <- 3
tol <- .3

test_that("Different input values for RA", {
  RAtest <- list(-1, 0, "0", NA)
  lapply(RAtest, function(x) {
    expect_error(do.call(x, racusum_beta_crit_sim, L0=L0, RQ=RQ, g0=g0, g1=g1, shape1=shape1, shape2=shape2, rs=maxS, RA=x))})
})

test_that("Different simulation algorithms, detecting deterioration", {
  skip_on_cran()
  skip_if(SKIP == TRUE, "skip this test now")

  ## RA=2
  m <- 1e3
  expect_equal(racusum_beta_crit_sim(L0=L0, RA=2, RQ=RQ, g0=g0, g1=g1, shape1=shape1, shape2=shape2, rs=maxS, verbose=TRUE, m=m), 2.5091, tolerance=tol)
  expect_equal(racusum_beta_crit_sim(L0=L0, RA=2, RQ=RQ, g0=g0, g1=g1, shape1=shape1, shape2=shape2, rs=maxS, verbose=FALSE, m=m), 2.5091, tolerance=tol)
})
