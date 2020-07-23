context("racusum_beta_FWT2")

L0 <- 500
g0 <- -3.6798
g1 <- 0.0768*71
shape1 <- 1
shape2 <- 3
QS <- 1
tol <- 10^-6

expect_equal(FWT2(w=-.5, QA=2, g0=g0, g1=g1, shape1=shape1, shape2=shape2, QS=QS), 0.002794664, tolerance=tol)
expect_equal(FWT2(w=.5,  QA=1/2, g0=g0, g1=g1, shape1=shape1, shape2=shape2, QS=QS), 0.9998782, tolerance=tol)
