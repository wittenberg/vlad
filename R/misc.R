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

#' @template search_delta

#' @author Philipp Wittenberg
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

#' @name VMASK3
#' @title Vmask3
#' @description Helper function to compute truncated symeterical/asymetrical vmask
#'
#' @param A ...
#' @param B ...
#' @param d1 Double. For the XYZ CUSUM Distance d from vertex of V-Mask. d=h/k
#' @param d2 Double. For the XYZ CUSUM Distance d from vertex of V-Mask. d=h/k
#' @param theta1 Double. Angle ...
#' @param theta2 Double. Angle ...
#' @param Sn ...
#' @param seg Logical. ...
#'
#' @return ...
#'
#' @importFrom tidyr gather

#' @author Philipp Wittenberg
#' @export
VMASK3 <- function(A, B, d1, d2, theta1, theta2, Sn, seg) {
  number <- n <- NULL
  theta1 <- (theta1*pi)/180
  theta2 <- (theta2*pi)/180
  ## helpers
  last.n <- B
  last.Sn <- Sn[last.n]
  upper <- function(i) last.Sn + tan(theta1) * (last.n - i + d1)
  lower <- function(i) last.Sn - tan(theta2) * (last.n - i + d2)
  if(seg==TRUE) { ## vertical segments only
    mask <- data.frame(cbind(n  =c(A, last.n),
                             up =c(upper(c(A, (last.n+d1)))),
                             low=c(lower(c(A, (last.n+d2))))),
                       number=rep(B, 2) )
    tidyr::gather(mask, key="group", value=value, c(-n, -number))
  } else { ## truncated V-mask
    mask <- data.frame(cbind(n  =c(A, last.n),
                             up =c(upper(c(A, (last.n)))),
                             low=c(lower(c(A, (last.n))))),
                       number=rep(B, 2) )
    tidyr::gather(mask, key="group", value=value, c(-n, -number))
  }
}

#' @name compute_vmask
#' @title Compute V-Masks arms, nose and alarm points
#' @description Function for plotting truncated symeterical/asymetrical vmask

#' @param z Numeric Vector. ...
#' @param d1 Double. For the XYZ CUSUM Distance d from vertex of V-Mask. d=h/k
#' @param d2 Double. For the XYZ CUSUM Distance d from vertex of V-Mask. d=h/k
#' @param theta1 Double. Angle ...
#' @param theta2 Double. Angle ...

#' @return ...

#' @importFrom dplyr bind_rows filter mutate select
#' @importFrom magrittr "%>%"

#' @author Philipp Wittenberg
#' @export
compute_vmask <- function(z, d1, d2, theta1, theta2) {
  Signal <- Clow <- Cup <- n <- k1 <- k2 <- h1 <- h2 <- NULL
  k1 <- tan(pi*theta1/180)
  k2 <- tan(pi*theta2/180)
  h1 <- k1*d1
  h2 <- k2*d2
  Sn <- cumsum(z)
  cv <- eocusum_scores(z=z, k1=k1, k2=k2, reset=TRUE, h1=h1, h2=h2)
  a  <- data.frame(cbind("n"=1:length(cv$s1), "Cup"=cv$s1, "Clow"=cv$s1l, "h2"=h2, "h1"=h1)) %>%
    dplyr::mutate(Signal=c(ifelse(Cup > h1, "Alarm U", ifelse(Clow < -h2, "Alarm L", "No Alarm")))) %>%
    dplyr::filter(Signal!="No Alarm") %>% dplyr::select(n)
  g <- nrow(a)
  ## arms of vmask
  DD1 <- dplyr::bind_rows(
    ## first Vmask
    VMASK3(1, a[1,], d1, d2, theta1, theta2, Sn, FALSE),
    ## Vmasks in between
    lapply(2:g, function(i) VMASK3(a[i-1,], a[i,], d1, d2, theta1, theta2, Sn, FALSE)),
    ## last Vmask
    VMASK3(a[g,], length(Sn), d1, d2, theta1, theta2, Sn, FALSE),
    .id="Masks"
  )
  ## vertical segments of truncated vmask
  DD2 <- dplyr::bind_rows(
    VMASK3(a[1,], a[1,], d1, d2, theta1, theta2, Sn, TRUE),
    lapply(2:g, function(i) VMASK3(a[i,], a[i,], d1, d2, theta1, theta2, Sn, TRUE)),
    VMASK3(length(Sn), length(Sn), d1, d2, theta1, theta2, Sn, TRUE),
    .id="Masks"
  )
  ## alarm points
  DD3 <- cbind(Masks=as.factor(1:g), n=a, value=Sn[a[, 1]])
  return(list("arms"=DD1, "nose"=DD2, "alarms"=DD3))
}
