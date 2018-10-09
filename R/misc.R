#' @name llr_cusum_scores
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
#' @template llr_cusum_scores
#'
#' @author Philipp Wittenberg
#' @export
llr_cusum_scores <- function(wt1, wt2, reset = FALSE, h1 = NULL, h2 = NULL) {
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
      s1[i+1] <- max(0, o1 + wt1[i])
      s1l[i+1] <- max(0, o2 + wt2[i])
    } else {
      o1 <- s1[i]
      o2 <- s1l[i]
      s1[i+1] <- max(0, o1 + wt1[i])
      s1l[i+1] <- max(0, o2 + wt2[i])
    }
  }
  return(list("s1" = s1[-1], "s1l" = s1l[-1]))
}
