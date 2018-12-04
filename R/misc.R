#' @name racusum_scores
#' @title Compute CUSUM scores based on the log-likelihood ratio statistic
#' @description Compute CUSUM scores based on the log-likelihood ratio statistic.
#'
#' @param wt1 Double. Log-likelihood ratio scores from function \code{\link{llr_score}}
#' for upper CUSUM.
#' @param wt2 Double. Log-likelihood ratio scores from function \code{\link{llr_score}}
#' for lower CUSUM.
#' @param reset Logical. If \code{FALSE} CUSUM statistic is not reset. If \code{TRUE} CUSUM
#' statistic is reset to \code{0} after a signal is issued.
#' @param h1 Double. Upper control limit of the CUSUM chart.
#' @param h2 Double. Lower control limit of the CUSUM chart.
#'
#' @return Returns a list with two components for the CUSUM scores.
#'
#' @template racusum_scores
#'
#' @author Philipp Wittenberg
#' @export
racusum_scores <- function(wt1, wt2, reset = FALSE, h1 = NULL, h2 = NULL) {
  n <- length(wt1)
  s1 <- rep(0, n)
  s1l <- s1
  for (i in 1:n) {
    if (reset == TRUE) {
      if (s1[i] > h1 | s1l[i] > h2) {
        o1 <- 0
        o2 <- 0
      } else {
        o1 <- s1[i]
        o2 <- s1l[i]
      }
      s1[i+1]  <- max(0, o1 + wt1[i])
      s1l[i+1] <- max(0, o2 + wt2[i])
    } else {
      o1 <- s1[i]
      o2 <- s1l[i]
      s1[i+1]  <- max(0, o1 + wt1[i])
      s1l[i+1] <- max(0, o2 + wt2[i])
    }
  }
  return(list("s1" = s1[-1], "s1l" = s1l[-1]))
}

#' @name eocusum_scores
#' @title Compute CUSUM scores based on E-O
#' @description Compute CUSUM scores based on E-O.
#'
#' @param z NumericVector. \code{E-O} values.
#' @param k1 Double. Reference value \code{k} for detecting improvement can be determined from function
#' \code{\link{optimal_k}}.
#' @param k2 Double. Reference value \code{k} for detecting deteroration can be determined from function
#' \code{\link{optimal_k}}.
#' @param reset Logical. If \code{FALSE} CUSUM statistic is not reset. If \code{TRUE} CUSUM
#' statistic is reset to \code{0} after a signal is issued.
#' @param h1 Double. Upper control limit of the CUSUM chart.
#' @param h2 Double. Lower control limit of the CUSUM chart.
#'
#' @return Returns a list with two components for the CUSUM scores.
#'
#' @template eocusum_scores
#'
#' @author Philipp Wittenberg
#' @export
eocusum_scores <- function(z, k1, k2, reset = FALSE, h1 = NULL, h2 = NULL) {
  n <- length(z)
  s1 <- rep(0, n)
  s1l <- s1
  for (i in 1:n) {
    if (reset == TRUE) {
      if (s1[i] > h1 | s1l[i] < -h2) {
        o1 <- 0
        o2 <- 0
      } else {
        o1 <- s1[i]
        o2 <- s1l[i]
      }
      s1[i+1]  <- max(0, o1 + z[i] - k1)
      s1l[i+1] <- min(0, o2 + z[i] + k2)
    } else {
      o1 <- s1[i]
      o2 <- s1l[i]
      s1[i+1]  <- max(0, o1 + z[i] - k1)
      s1l[i+1] <- min(0, o2 + z[i] + k2)
    }
  }
  return(list("s1" = s1[-1], "s1l" = s1l[-1]))
}

#' @name search_delta
#' @title Search Box-Cox transformation parameter
#' @description Search Box-Cox transformation parameter.
#'
#' @param s Integer vector. Parsonnet Score values within a range of \code{0} to \code{100}
#' representing the preoperative patient risk.
#' @param y Double. Binary (\code{0/1}) outcome values of each operation.
#' @param type Character. If \code{type = "ML"} Maximum Likelihood used to search the Box-Cox
#'  transformation parameter, \code{type = "Pearson"} uses a Pearson measure.
#' @param dmin Double. Minimum value for the grid search.
#' @param dmax Double. Maximum value for the grid search.
#' @return Returns a single value for the Box-Cox transformation parameter.

#' @author Philipp Wittenberg
#' @examples
#' \dontrun{
#' ## load data
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## preprocess data to 30 day mortality and subset data to
#' ## phase I (In-control) and phase II (monitoring)
#' SALL <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II")))
#'
#' ## subset phase I (In-control)
#' SI <- filter(SALL, phase == "I") %>% select(s, y)
#'
#' ## search delta
#' dML <- search_delta(SI$s, SI$y, type = "ML")
#' dQQ <- search_delta(SI$s, SI$y, type = "Pearson")
#'
#' ## show Log-likelihood (ell()) and Pearson measure (QQ()) for each delta
#' delta <- c(-2, -1, 0, dML, dQQ, 0.5, 1, 2)
#' r <- sapply(delta, function(i) rbind(i, ell(SI$s, SI$y, i), QQ(SI$s, SI$y, i)))
#' rownames(r) <- c("d", "l", "S")
#' t(r)
#' data.frame(t(r)) %>% filter(l == max(l) | S == min(S))
#' }
#' @export
search_delta <- function(s, y, type = "ML", dmin = -2, dmax = 2) {
  switch(type,
         ML = {
           fmax  <- Vectorize(function(delta) - as.numeric( stats::logLik( stats::glm(y ~ trafo(delta, s), family = stats::binomial(link = "logit")) ) ) )
           delta <- stats::optimize(fmax, c(dmin, dmax), tol = 1e-9)$minimum
           delta
         },
         Pearson={
          fmin  <- Vectorize(function(delta) QQ(s, y, delta))
          delta <- stats::optimize(fmin, c(dmin, dmax, tol = 1e-9))$minimum
          delta
        }
  )
  delta
}

#' @name ell
#' @title Estimated log-likelihood.
#' @description Estimated log-likelihood.
#'
#' @param s Integer vector. Parsonnet Score values within a range of \code{0} to \code{100}
#' representing the preoperative patient risk.
#' @param y Double. Binary (\code{0/1}) outcome values of each operation.
#' @param delta Double. Box-Cox transformation parameter.
#' @return Returns a single value which is estimated log-likelihood.

#' @author Philipp Wittenberg
#' @examples
#' \dontrun{
#' ## load data
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## preprocess data to 30 day mortality and subset data to
#' ## phase I (In-control) and phase II (monitoring)
#' SALL <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II")))
#'
#' ## subset phase I (In-control)
#' SI <- filter(SALL, phase == "I") %>% select(s, y)
#'
#' dML <- search_delta(SI$s, SI$y, type = "ML")
#' ell(SI$s, SI$y, dML)
#' }
#' @export
ell <- function(s, y, delta) {
  stats::logLik(stats::glm(y ~ trafo(delta, s), family = stats::binomial(link = "logit")))
}

#' @name QQ
#' @title Pearson measure
#' @description Pearson measure.
#'
#' @param s Integer vector. Parsonnet Score values within a range of \code{0} to \code{100}
#' representing the preoperative patient risk.
#' @param y Numeric Vector. Binary (\code{0/1}) outcome values of each operation.
#' @param delta Double. Box-Cox transformation parameter.
#' @return Returns a single value.

#' @author Philipp Wittenberg
#' @examples
#' \dontrun{
#' ## load data
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## preprocess data to 30 day mortality and subset data to
#' ## phase I (In-control) and phase II (monitoring)
#' SALL <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II")))
#'
#' ## subset phase I (In-control)
#' SI <- filter(SALL, phase == "I") %>% select(s, y)
#'
#' dQQ <- search_delta(SI$s, SI$y, type = "Pearson")
#' QQ(SI$s, SI$y, dQQ)
#' }
#' @export
QQ <- function(s, y, delta) {
  GLM <- stats::glm(y ~ trafo(delta, s), family = stats::binomial(link = "logit"))
  alpha <- stats::coef(GLM)[1]
  beta  <- stats::coef(GLM)[2]
  pi.of.s <- function(s) 1 / (1 + exp(-alpha - beta * trafo(delta, s)))
  s_ <- sort(unique(s))
  pi.hat <- pi.of.s(s_)
  ybar <- tapply(y, s, mean)
  nn <- tapply(y, s, length)
  QQ <- sum(nn * (pi.hat - ybar)^2 / pi.hat / (1 - pi.hat))
  QQ
}

#' @name trafo
#' @title Box-Cox transformation of data.
#' @description Box-Cox transformation of data.
#'
#' @param delta Numeric. Box-Cox transformation parameter.
#' @param x Numceric Vector. Parsonnet Score values within a range of 0 to 100 representing the
#' preoperative patient risk.
#' @return Returns a transformed Numeric vector.

#' @author Philipp Wittenberg
#' @export
trafo <- Vectorize(function(delta, x) ifelse( abs(delta)<1e-9, log(1+x), ((1+x)^delta-1)/delta ) )
