
## redirect graphical output to pdf file
pdf("Supplemental-material-SIM-18-0571.pdf", width=8, height=6, pointsize=8)

require(vlad)

# data transferred from the R package "scpadjust" version 1.1 (2015-11-20): Axel Gandy and Jan Terje Kvaloy
#
#require(spcadjust)
#data("cardiacsurgery", package = "spcadjust")
#y <- as.numeric( cardiacsurgery$time <= 30  )
#cardiacsurgery <- cbind(cardiacsurgery, y)
#
#N0 <- table( cardiacsurgery$Parsonnet )
#S0 <- as.numeric( names(N0) )
#N0 <- as.numeric(N0)
#X0 <- as.numeric( tapply(cardiacsurgery$y, cardiacsurgery$Parsonnet, sum) )
#DD0 <- data.frame(S0, N0, X0)
#
#cardiacsurgery <- cardiacsurgery[1:1766,]
#
#N1 <- table( cardiacsurgery$Parsonnet )
#S1 <- as.numeric( names(N1) )
#N1 <- as.numeric(N1)
#X1 <- as.numeric( tapply(cardiacsurgery$y, cardiacsurgery$Parsonnet, sum) )
#DD1 <- data.frame(S1, N1, X1)
#

# the same data, stored explicitly
#
## Gandy and Kvaloy 2015, complete data
DD0 <- data.frame(S0=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 65, 67, 69, 71),
                  N0=c(850, 147, 330, 560, 222, 399, 258, 340, 186, 182, 238, 151, 196, 129, 114, 139, 95, 118, 83, 144, 80, 57, 59, 49, 66, 43, 36, 33, 27, 31, 25, 13, 21, 12, 10, 11, 5, 6, 7, 5, 12, 4, 7, 3, 5, 4, 7, 6, 2, 1, 9, 6, 5, 7, 5, 4, 5, 5, 5, 3, 3, 1, 4, 1, 2, 1, 1),
                  X0=c(7, 3, 6, 9, 6, 13, 10, 16, 10, 11, 14, 10, 12, 19, 9, 11, 5, 12, 11, 20, 5, 4, 14, 7, 12, 5, 6, 7, 2, 3, 7, 1, 6, 2, 3, 2, 1, 3, 1, 3, 3, 1, 2, 0, 2, 2, 3, 0, 1, 1, 5, 3, 2, 5, 5, 1, 1, 3, 3, 3, 2, 1, 2, 1, 1, 0, 0))

## Gandy and Kvaloy 2015, first two years (phase I)
DD1 <- data.frame(S1=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 48, 49, 50, 51, 52, 53, 54, 57, 58, 60, 62, 67, 69),
                  N1=c(312, 53, 131, 166, 88, 128, 74, 106, 54, 65, 68, 51, 57, 23, 32, 27, 21, 24, 25, 51, 27, 15, 11, 9, 15, 10, 15, 12, 8, 11, 4, 5, 10, 3, 3, 7, 2, 1, 1, 2, 3, 4, 3, 1, 1, 2, 3, 1, 1, 2, 3, 2, 2, 1, 1, 3, 2, 2, 1, 1),
                  X1=c(1, 0, 0, 3, 4, 3, 2, 2, 5, 4, 5, 2, 7, 5, 3, 3, 2, 1, 4, 8, 0, 2, 0, 2, 5, 0, 3, 5, 1, 0, 1, 0, 4, 0, 0, 2, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 2, 1, 0, 2, 1, 1, 1, 0))

attach(DD0)
attach(DD1)

ps <- 0:71

### Figure 1 ###

plot(N0, type="h", xlab="Parsonnet score s", ylab="Frequency", main="Figure 1 (A)")


# original model, all data
GLM0a <- glm(cbind(X0,N0-X0) ~ S0, family="binomial")
P0a <- predict(GLM0a, newdata=data.frame(S0=ps), type="response")

# log score model, all data
lS0 <- log(S0+1)
GLM0b <- glm(cbind(X0,N0-X0) ~ lS0, family="binomial")
P0b <- predict(GLM0b, newdata=data.frame(lS0=log(ps+1)), type="response")

# original model, first two years (phase 1)
GLM1 <- glm(cbind(X1,N1-X1) ~ S1, family="binomial")
P1 <- predict(GLM1, newdata=data.frame(S1=ps), type="response")

H0 <- X0 / N0

plot(S0, H0, pch=19, xlab="Parsonnet score s", ylab="Estimated probability of death", main="Figure 1 (B)")

lines(ps, P0a, lwd=2, lty=1, col="red")
lines(ps, P1, lwd=2, lty=2, col="red")
lines(ps, P0b, lwd=2, lty=4, col="blue")

legend("topleft", c("logit, 2 years", "logit, all data", "logit(log), all data"), lty=c(2, 1, 4), col=c("red", "red", "blue"), lwd=3, bg="white", seg.len=4)


### Figure 2 ###

# original model, all patients
plot(ps, P0a, type="l", lwd=2, xlab="Parsonnet score s", ylab="Estimated probability of death", main="Figure 2 (A)", ylim=c(0,1))

LTY <- COL <- 2
for ( smax in c(64, 61, 58, 56, 54, 52) ) { # remove patients not better than smax
  s0 <- S0[S0 < smax]
  n0 <- N0[S0 < smax]
  x0 <- X0[S0 < smax]
  GLM0a2 <- glm(cbind(x0,n0-x0) ~ s0, family="binomial")
  P0a2 <- predict(GLM0a2, newdata=data.frame(s0=ps), type="response")
  lines(ps, P0a2, lty=LTY, col=COL, lwd=2)
  LTY <- LTY + 1
  COL <- COL + 1
}

legend("topleft", c("all", "<64", "<61", "<58", "<56", "<54", "<52"), lwd=3, lty=1:7, col=1:7, bg="white", seg.len=4)


# logit(log) model, all patients
plot(ps, P0b, type="l", lwd=2, xlab="Parsonnet score s", ylab="Estimated probability of death", main="Figure 2 (B)", ylim=c(0,1))

LTY <- COL <- 2
for ( smax in c(64, 61, 58, 56, 54, 52) ) { # remove patients not better than smax
  s0 <- S0[S0 < smax]
  n0 <- N0[S0 < smax]
  x0 <- X0[S0 < smax]
  ls0 <- log(s0+1)
  GLM0b2 <- glm(cbind(x0,n0-x0) ~ ls0, family="binomial")
  P0b2 <- predict(GLM0b2, newdata=data.frame(ls0=log(ps+1)), type="response")
  lines(ps, P0b2, lty=LTY, col=COL, lwd=2)
  LTY <- LTY + 1
  COL <- COL + 1
}

legend("topleft", c("all", "<64", "<61", "<58", "<56", "<54", "<52"),  lwd=3, lty=1:7, col=1:7, bg="white", seg.len=4)


### Figure 3 ###

D <- c(-2, -1, 0, .249, .5, 1, 2)

# response profile, original model
plot(ps, P1, type="l", lwd=2, xlab="Parsonnet score s", ylab="Estimated probability of death", main="Figure 3 (A)", ylim=c(0,1), lty=6, col=6)
abline(v=0, h=0, lty=2, col="grey")

LTY <- COL <- 1
for ( delta in D ) { # loop over various Box-Cox delta
  ds1 <- trafo(delta, S1)
  GLM1 <- glm(cbind(X1,N1-X1) ~ ds1, family="binomial")
  P1  <- predict(GLM1, newdata=data.frame(ds1=trafo(delta, ps)), type="response")
  lines(ps, P1, lty=LTY, col=COL, lwd=2)
  LTY <- LTY + 1
  COL <- COL + 1
}

legend("topleft", as.character(D), lwd=3, lty=1:7, col=1:7, bg="white", seg.len=4, title=expression(delta))


# W for surviving patient, detect deterioriation
QA <- 2
plot(ps, -log(1-P1+QA*P1), type="l", lwd=2, xlab="Parsonnet score s", ylab="Estimated probability of death", main="Figure 3 (B)", ylim=c(-0.8,0), lty=6, col=6)

abline(v=0, h=0, lty=2, col="grey")

LTY <- COL <- 1
for ( delta in D ) { # loop over various Box-Cox delta
  ds1 <- trafo(delta, S1)
  GLM1 <- glm(cbind(X1,N1-X1) ~ ds1, family="binomial")
  P1  <- predict(GLM1, newdata=data.frame(ds1=trafo(delta, ps)), type="response")
  lines(ps, -log(1-P1+QA*P1), lty=LTY, col=COL, lwd=2)
  LTY <- LTY + 1
  COL <- COL + 1
}

legend("bottomleft", as.character(D), lwd=3, lty=1:7, col=1:7, bg="white", seg.len=4, title=expression(delta))


# W for surviving patient, detect improvement
QA <- 1/2
plot(ps, -log(1-P1+QA*P1), type="l", lwd=2, xlab="Parsonnet score s", ylab="Estimated probability of death", main="Figure 3 (C)", ylim=c(0,0.8), lty=6, col=6)

abline(v=0, h=0, lty=2, col="grey")

LTY <- COL <- 1
for ( delta in D ) { # loop over various Box-Cox delta
  ds1 <- trafo(delta, S1)
  GLM1 <- glm(cbind(X1,N1-X1) ~ ds1, family="binomial")
  P1  <- predict(GLM1, newdata=data.frame(ds1=trafo(delta, ps)), type="response")
  lines(ps, -log(1-P1+QA*P1), lty=LTY, col=COL, lwd=2)
  LTY <- LTY + 1
  COL <- COL + 1
}

legend("topleft", as.character(D), lwd=3, lty=1:7, col=1:7, bg="white", seg.len=4, title=expression(delta))


### Figure 4 ###

# calculate log-likelihood, correction applied for grouped data
lol <- Vectorize(function(delta) {
  ds1 <- trafo(delta, S1)
  GLM <- glm(cbind(X1,N1-X1) ~ ds1, family="binomial")
  korr <- sum( log(choose(N1, X1)) )
  logLik(GLM) - korr
})

curve(lol, -2, 2, xlab=expression(delta), ylab="log-likelihood", main="Figure 4")

OPT <- optimize(lol, c(-2, 2), maximum=TRUE)

abline(v=OPT$maximum, h=OPT$objective, lty=2, col="grey")

points(OPT$maximum, OPT$objective, pch=19, col="red")


### Table 1 ###

# calculate Pearson measure
chi <- Vectorize(function(delta) {
  ds1 <- trafo(delta, S1)
  GLM <- glm(cbind(X1,N1-X1) ~ ds1, family="binomial")
  P   <- predict(GLM, newdata=data.frame(ds1=ds1), type="response")
  S   <- sum( (P - X1/N1)^2 / ( P*(1-P) ) * N1 )
  S
})

D <- c(-2, -1, -0.5, 0, 0.063, .249, .256, .5, 1, 2)

print( cbind(D, round(lol(D), digits=3), round(chi(D), digits=3)) )


### Figure 5 ###

# only calculation of the CUSUM trhesholds

sca <- 3000
L0 <- 3800

# delta = 1 alias original model

GLM <- glm(cbind(X1,N1-X1) ~ S1, family="binomial")
P <- predict(GLM, newdata=data.frame(S1=S1), type="response")
pmix <- data.frame(h=N1/sum(N1), p1=P, p2=P)

h1 <- racusum_crit_mc(pmix=pmix, L0=L0, scaling=sca, RQ=1, RA=2)
h2 <- racusum_crit_mc(pmix=pmix, L0=L0, scaling=sca, RQ=1, RA=1/2)

cat(paste("\n\nThresholds for original model: h.upper =", h1, " and h.lower =", h2, "\n\n"))


# delta = 0.063 alias 'optimized' model (different to the optimal delta value for the original data)

delta <- 0.063
dS1 <- trafo(delta, S1)
GLM <- glm(cbind(X1,N1-X1) ~ dS1, family="binomial")
P <- predict(GLM, newdata=data.frame(dS1=dS1), type="response")
pmix <- data.frame(h=N1/sum(N1), p1=P, p2=P)

h3 <- racusum_crit_mc(pmix=pmix, L0=L0, scaling=sca, RQ=1, RA=2)
h4 <- racusum_crit_mc(pmix=pmix, L0=L0, scaling=sca, RQ=1, RA=1/2)

cat(paste("\n\nThresholds for optimized model (delta=0.063): h.upper =", h3, " and h.lower =", h4, "\n\n"))


### Figure 6 ###

sca <- 3000
L0 <- 740

# original model, delta = 1
oGLM <- glm(cbind(X1,N1-X1) ~ S1, family="binomial")

F6 <- function(smin, QA=2, LTY=1, SF="(A)") { # remove healthy patients with score < smin
  n1 <- N1[S1 >= smin]
  s1 <- S1[S1 >= smin]
  P <- predict(oGLM, newdata=data.frame(S1=s1), type="response") # assumed model (delta = 1)
  pmix <- data.frame(h=n1/sum(n1), p1=P, p2=P)
  h <- racusum_crit_mc(pmix=pmix, L0=L0, scaling=sca, RQ=1, RA=QA)

  ARL0 <- Vectorize(function(delta) {
    dS1 <- trafo(delta, S1)
    GLM <- glm(cbind(X1,N1-X1) ~ dS1, family="binomial") # true model, mostly different to original one
    PP  <- predict(GLM, newdata=data.frame(dS1=trafo(delta, s1)), type="response")
    pmix <- data.frame(h=n1/sum(n1), p1=P, p2=PP)
    L0 <- racusum_arl_mc(pmix=pmix, RA=QA, RQ=1, h=h, scaling=sca)
    L0
  })

  if ( smin==0 ) {
    curve(ARL0, -2, 2, xlab=expression(delta), ylab=expression(ARL[0]), ylim=c(300,1100), n=51, lwd=2, main=paste("Figure 6", SF))
    abline(v=c(.063, 1), h=L0, col="grey", lty=2)
  } else {
    curve(ARL0, -2, 2, add=TRUE, n=51, lwd=2, lty=LTY)
  }
}


F6(0)
F6(2, LTY=2)
F6(6, LTY=3)

legend("topleft", c(expression(D[100]), expression(D[79]), expression(D[50])), lty=1:3, lwd=3, seg.len=4, bg="white")


F6(0, QA=0.5, SF="(B)")
F6(2, QA=0.5, LTY=2)
F6(6, QA=0.5, LTY=3)

legend("topleft", c(expression(D[100]), expression(D[79]), expression(D[50])), lty=1:3, lwd=3, seg.len=4, bg="white")


### Table 2 ###

n1 <- N1[S1 >= 2]
s1 <- S1[S1 >= 2]
P <- predict(oGLM, newdata=data.frame(S1=s1), type="response") # assumed model (delta = 1)
pmix <- data.frame(h=n1/sum(n1), p1=P, p2=P)
h1 <- racusum_crit_mc(pmix=pmix, L0=L0, scaling=sca, RQ=1, RA=2)
h2 <- racusum_crit_mc(pmix=pmix, L0=L0, scaling=sca, RQ=1, RA=1/2)

ARL0 <- Vectorize(function(delta, QA, h) {
  dS1 <- trafo(delta, S1)
  GLM <- glm(cbind(X1,N1-X1) ~ dS1, family="binomial") # true model, mostly different to original one
  PP  <- predict(GLM, newdata=data.frame(dS1=trafo(delta, s1)), type="response")
  pmix <- data.frame(h=n1/sum(n1), p1=P, p2=PP)
  L0 <- racusum_arl_mc(pmix=pmix, RA=QA, RQ=1, h=h, scaling=sca)
  L0
}, "delta")

D <- c(-2, -1, 0, 0.063, .25, .5, 1, 2)

print( cbind(D, round(ARL0(D, 2, h1), digits=2), round(ARL0(D, 1/2, h2), digits=2)) )


### Figure A1 ###

# Steiner et al. (2000) setup (detect deterioriation)
hS <- 4.5
QA <- 2

GLM <- glm(cbind(X1,N1-X1) ~ S1, family="binomial") # original model, first 2 years
P <- predict(GLM, newdata=data.frame(S1=S1), type="response")
pmix <- data.frame(h=N1/sum(N1), p1=P, p2=P)

# calculate 'asymptotic result (gamma = 10,000)
LLlarge <- racusum_arl_mc(hS, pmix, QA, 1, scaling=1e4, rounding="p", method="Toep")

# loop over many scaling values
gg <- seq(50, 1000, by=1)
LLs <- LLp <- rep(NA, length(gg))
for ( i in 1:length(gg) ) {
  LLs[i] <- racusum_arl_mc(hS, pmix, QA, 1, scaling=gg[i], rounding="s", method="Toep")
  LLp[i] <- racusum_arl_mc(hS, pmix, QA, 1, scaling=gg[i], rounding="p", method="Toep")
}


plot(gg, LLs, type="l", xlab=expression(gamma), ylab="ARL approximation", ylim=LLlarge+c(-1,1)*3e3, main="Figure A1 (A)")

lines(gg, LLp, col="blue")

legend("topright", c("Steiner", "rounding pairs"), lty=1, lwd=3, col=c("black", "blue"), bg="white", seg.len=4)


plot(gg[1+(0:95)*10], LLp[1+(0:95)*10], type="l", xlab=expression(gamma), ylab="ARL approximation", ylim=LLlarge+c(-1,1)*50, col="blue", main="Figure A1 (B)")

abline(h=LLlarge, col="red", lty=2, lwd=2)

legend("topright", c("rounding pairs", expression(gamma == 10^4)), lty=1:2, lwd=3, col=c("blue", "red"), bg="white", seg.len=4)


dev.off()
