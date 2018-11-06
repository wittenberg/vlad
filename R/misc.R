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

#' @name rockettails
#' @title Rockettails for VLAD
#' @description Rockettails for VLAD.
#'
#' @param alpha double.
#' @param s double. Parsonnet score.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model. For more information see details.
#'
#' @return ....
#'
#' @keywords keyword1 keyword2
#' @author Philipp Wittenberg
#' @examples
#' \dontrun{
#' library("dplyr")
#' library("tidyr")
#' library("ggplot2")
#' data("cardiacsurgery", package = "spcadjust")
#' SALL <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II")))
#'
#' SALLI <- filter(SALL, phase == "I")
#' coeff1 <- round(coef(glm(y ~ s, data = SALLI, family = "binomial")), 3)
#'
#' S2 <- filter(SALL, surgeon == 2)
#' S2I <- filter(S2, phase == "I") %>% select(s, y)
#' S2II <- filter(S2, phase == "II") %>% select(s, y)
#'
#' rt <- rockettails(alpha = 0.05, s = S2II[, "s"] , coeff = coeff1)
#' EO <- sapply(1:nrow(S2II), function(i) calceo(df = S2II[i, c("s", "y")], coeff = coeff1,
#'              yemp = TRUE))
#' df1 <- data.frame(cbind("n" = 1:nrow(S2II), "cEO" = cumsum(EO), rt)) %>%
#'   gather(key = variable, value = value, c(-n))
#'
#' ggplot(df1, aes(x = n, y = value, group = variable)) +
#'   geom_line() +
#'   geom_point(data = filter(df1, variable == "cEO"), col = "red") +
#'   geom_hline(yintercept = 0, linetype = "dashed", col = "darkgreen") +
#'   theme_bw() + labs(x = "Patient number n", y = "CUSUM E-O") +
#'   geom_hline(yintercept=0, linetype="dashed", col="darkgreen") +
#'   scale_y_continuous(sec.axis=dup_axis(name=NULL, labels=NULL)) +
#'   scale_x_continuous(sec.axis=dup_axis(name=NULL, labels=NULL))
#' }
#' @export
rockettails <- function(alpha, s, coeff){
  y <- sapply(1:length(s), function(i) gettherisk(s[i], coeff = coeff))
  r <- sapply(1:length(s), function(i) y[i]*(1-y[i]))
  cbind("upper tail" = sqrt(cumsum(r))*stats::qnorm(1-alpha/2, lower.tail = TRUE),
        "lower tail" = sqrt(cumsum(r))*stats::qnorm(1-alpha/2, lower.tail = FALSE)
  )
}

#' @name plot_racusum
#' @title Plot a two-sided risk-adjusted CUSUM chart
#' @description Plot a two-sided risk-adjusted CUSUM chart.
#'
#' @param data Dataframe. TODO
#' @param coeff Numeric Vector. TODO
#' @param signal Logical. TODO
#' @inheritParams racusum_scores
#'
#' @return Returns a ggplot object for a two-sided risk-adjusted CUSUM chart.
#'
#' @importFrom tidyr gather
#' @import dplyr
#' @import ggplot2
#'
#' @template plot_racusum
#'
#' @author Philipp Wittenberg
#' @export
plot_racusum <- function(data, coeff, h1, h2, reset = FALSE, signal = FALSE) {
  ## calculate CUSUM weights
  wt1 <- sapply(1:nrow(data), function(i) llr_score(data[i, c("s", "y")], coeff = coeff, RA = 2))
  wt2 <- sapply(1:nrow(data), function(i) llr_score(data[i, c("s", "y")], coeff = coeff, RA = 1/2))
  ## CUSUM statistics with/without reset
  cv <- racusum_scores(wt1 = wt1, wt2 = wt2, reset = reset, h1 = h1, h2 = h2)
  ## CUSUM values and limits; determine signals
  dm1 <- data.frame(cbind("n"    = 1:length(wt1),
                          "Cup"  = cv$s1,
                          "Clow" = -cv$s1l,
                          "UCL"   = h1,
                          "LCL"   = -h2)) %>%
    mutate("n" = row_number()) %>%
    gather("CUSUM", value, -n) %>%
    mutate(signal = ifelse(CUSUM == "Clow" & value < -h2, "AlarmL",
                           ifelse(CUSUM == "Cup" & value > h1, "AlarmU", "No Alarm")))

  dm1$CUSUM <- factor(dm1$CUSUM, levels = c("UCL", "Cup", "Clow", "LCL"))
  labl <- list("UCL", expression(C[up]), expression(C[low]), "LCL")

  ## plot
  p <- ggplot(dm1, aes(x = n, y = value, colour = CUSUM, group = CUSUM)) +
    geom_hline(yintercept = 0, colour = "darkgreen", linetype = "dashed") +
    geom_line(size = 0.5) +
    labs(x = "Patient number n", y = "CUSUM values") +
    theme_classic() +
    scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
    scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
    scale_color_manual(values = c("red", "orange", "blue", "red"), labels = labl,
                       guide_legend(title = ""))
  ## add alarm signals as points
  if (signal == TRUE) {
    p + geom_point(data = filter(dm1, signal == "AlarmU" | signal == "AlarmL"),
                     colour = "red")
  } else p
}
