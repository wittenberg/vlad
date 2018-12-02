setwd("~/clones/vlad//data-raw")
######### Load data #############################
load("cardiacsurgery.rda")
y <- as.numeric( cardiacsurgery$time <= 30  )
cardiacsurgery <- cbind(cardiacsurgery, y)

N2 <- table( cardiacsurgery$Parsonnet )
S2 <- as.numeric( names(N2) )
N2 <- as.numeric(N2)
X2 <- as.numeric( tapply(cardiacsurgery$y, cardiacsurgery$Parsonnet, sum) )

GandyKvaloy_all <- data.frame(S0=S2, N0=N2, X0=X2)
usethis::use_data(GandyKvaloy_all, overwrite = TRUE)

cardiacsurgery <- subset(cardiacsurgery, c(date <= 730))

N3 <- table( cardiacsurgery$Parsonnet )
S3 <- as.numeric( names(N3) )
N3 <- as.numeric(N3)
X3 <- as.numeric( tapply(cardiacsurgery$y, cardiacsurgery$Parsonnet, sum) )

GandyKvaloy_P1 <- data.frame(S1=S3, N1=N3, X1=X3)
usethis::use_data(GandyKvaloy_P1, overwrite = TRUE)
