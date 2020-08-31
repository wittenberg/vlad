## set up random number generation
set.seed(seed=1234)

## risk model, parameter vector
b <- c(-3.6798, .0768)

## Beta-binomial distribution parameters
n <- 83
ab <- c(shape1=.59, shape2=4.12)

## sample size
N <- 2500

## sample risk scores
S <- sapply(1:N, function(i, n, ab) rbinom(1, n, rbeta(i, ab[1], ab[2])), n=n, ab=ab)

## generate binary outcome
outgen <- function(b, s, RQ) {
  logitp <- b[1] + s * b[2]
  pt     <- exp(logitp)/(1 + exp(logitp))
  xstar  <- (RQ*pt) / (1-pt+RQ*pt)
  y <- as.numeric( runif(1) < xstar )
  c(s, y)
}
DF <- sapply(S[1:N], function(i) outgen(b, i, RQ=1))
surgery <- data.frame(s=DF[1,], y=DF[2,])

## save data
usethis::use_data(surgery, overwrite=FALSE)
