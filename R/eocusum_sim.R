#' @name optimal.k
#' @title Compute optimal k
#' @description Compute optimal k.
#'
#' @param QA double. Defines the performance of a surgeon with the odds ratio ratio of death Q.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model. For more information see details.
#' @param parsonnetscores NumericVector. Vector of Parsonnet Scores.
#'
#' @return output A description of the object the function outputs
#'
#' @details Formula deterioration:  \deqn{ k{det} = \frac{Q{A} - 1 - log(Q{A})}{log(Q{A})}\bar{p} , Q{A} > 1    }
#'          Formula improvement:    \deqn{ k{imp} = \frac{1 - Q{A} + log(Q{A})}{log(Q{A})}\bar{p} , Q{A} < 1    }
#'
#' @author Philipp Wittenberg
#' @examples
#' \dontrun{
#' library("vlad"); library("spcadjust")
#' data("cardiacsurgery")
#' cardiacsurgery <- dplyr::mutate(cardiacsurgery, phase=factor(ifelse(date < 2*365, "I", "II")))
#' S2I <- subset(cardiacsurgery, c(surgeon==2 & phase=="I"), c("Parsonnet", "status"))
#' coeff <- coef(glm(status ~ Parsonnet, data=S2I, family="binomial"))
#' kopt <- optimal.k(QA=2, parsonnetscores=S2I$Parsonnet, coeff=coeff)
#' kopt ## (Deterioration)
## k_opt = 0.04059649
#' kopt <- optimal.k(QA=1/2, parsonnetscores=S2I$Parsonnet, coeff=coeff)
#' kopt ##(Improvement)
#' ## k_opt = 0.02555328
#' ### k = kopt
#' QA <- 1/2
#' # manually find optimal k
#' pbar <- mean(sapply(df1[, 1], gettherisk, coef=coeff))
#' kopt <- pbar * ( 1 - QA + log(QA) ) / log(QA)
#' all.equal(kopt, optimal.k(QA=1/2, parsonnetscores=S2I$Parsonnet, coeff=coeff) )
#' }
#' @export
optimal.k <- function(QA, parsonnetscores, coeff) {
  .optimal_k(
    as.numeric(QA),
    as.vector(parsonnetscores),
    as.vector(coeff)
  )
}

#' @name gettherisk
#' @title Compute Risk of death
#' @description Compute Risk of death.
#'
#' @param parsonnetscore int. Parsonnet Score.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model. For more information see details.
#'
#' @return output A description of the object the function outputs
#'
#' @template gettherisk
#'
#' @author Philipp Wittenberg
#'
#' @export
gettherisk <- function(parsonnetscore, coeff) {
  .gettherisk(
    as.integer(parsonnetscore),
    as.vector(coeff)
  )
}

#' @name calceo
#' @title Compute Expected minus Observed value
#' @description Compute Expected minus Observed value.
#'
#' @param df DataFrame. First column Parsonnet Score and second column outcome of each operation.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model. For more information see details.
#' @param yemp boolean. If TRUE use observed outcome value, if FALSE use estimated binary logistc
#'  regression model.
#'
#' @return output A description of the object the function outputs
#'
#' @template calceo
#'
#' @keywords keyword1 keyword2
#' @family VLAD functions
#' @author Philipp Wittenberg
#' @export
calceo <- function(df, coeff, yemp = TRUE) {
  .calceo(
    as.data.frame(df),
    as.vector(coeff),
    as.logical(yemp)
  )
}

#' @name eocusum.arl.sim
#' @title Compute ARLs of EO-CUSUM control charts using simulation
#' @description Compute ARLs of EO-CUSUM control charts using simulation.
#'
#' @param r int. Number of of simulation runs.
#' @param k double. Reference value of the CUSUM control chart.
#' @param h double. Decision interval (alarm limit, threshold) of the CUSUM control chart.
#' @inheritParams calceo
#' @param side int. If side = 1, calculate ARL for the upper arm of the V-mask. If side = 2,
#'  calulate the lower arm of the V-mask.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @details Describe the algorithm for calulating in control run length ...
#'
#' @template eocusum.arl.sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum.arl.sim <- function(r, k, h, df, coeff, yemp = TRUE, side = 1) {
  .eocusum_arl_sim(
    as.integer(r),
    as.numeric(k),
    as.numeric(h),
    as.data.frame(df),
    as.vector(coeff),
    as.logical(yemp),
    as.integer(side)
  )
}

#' @name eocusum.arloc.sim
#' @title Compute Out of Control ARLs of EO-CUSUM control charts using simulation
#' @description Compute Out of Control ARLs of EO-CUSUM control charts using simulation.
#'
#' @inheritParams eocusum.arl.sim
#' @param coeff2 NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model of a resampled dataset. For more information see
#'  details.
#' @param QS double. Defines the performance of a surgeon with the odds ratio ratio of death Q.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @details Describe the algorithm for calulating in control run length ...
#'
#' @template eocusum.arloc.sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum.arloc.sim <- function(r, k, h, df, coeff, coeff2, QS = 1, side = 1) {
  .eocusum_arloc_sim(
    as.integer(r),
    as.numeric(k),
    as.numeric(h),
    as.data.frame(df),
    as.vector(coeff),
    as.vector(coeff2),
    as.numeric(QS),
    as.integer(side)
  )
}

#' @name eocusum.adoc.sim
#' @title Compute steady-state ARLs of EO-CUSUM control charts using simulation
#' @description Compute steady-state ARLs of EO-CUSUM control charts using simulation.
#'
#' @inheritParams eocusum.arloc.sim
#' @param m Integer. Simulated in-control observations.
#' @param type character. Default argument is "cond" for computation of conditional steady-state.
#' Other option is the cyclical steady-state "cycl".
#'
#' @return Returns a single value which is the Run Length.
#'
#' @details Describe the resampling algorithm for calulating out of control run length.
#'
#' @template eocusum.adoc.sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum.adoc.sim <- function(r, k, h, df, coeff, coeff2, QS = 1, side = "low", type = "cond", m = 50) {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("number of simulation runs r must a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("reference value k must a positive numeric value")}
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("control limit h must a positive numeric value")}
  df <- as.data.frame(df)
  if (class(df) != "data.frame") {stop("provide a dataframe for argument \"df\"")}
  else if (ncol(df) != 2) {stop("provide a dataframe with two columns for argument \"df\"")}
  else if (sapply(df, class)[1] != "integer") {stop("first column of dataframe must be of type integer")}
  else if (sapply(df, class)[2] != "numeric") {stop("second column of dataframe must be of type numeric")}
  coeff <- as.vector(coeff)
  if (is.na(coeff)  || length(coeff)  != 2) {stop("model coefficients \"coeff\"  must a numeric vector with two elements")}
  if (is.na(coeff2) || length(coeff2) != 2) {stop("model coefficients \"coeff2\" must a numeric vector with two elements")}
  iside <- switch(as.character(side), low = 1, up = 2)
  if (is.null(iside)) {
    warning("no valid input, using side=low (deterioration) as default")
    iside <- 1
  }
  QS <- as.numeric(QS)
  if (is.na(QS) || QS < 0) {stop("QS must a positive numeric value")}
  else if (QS < 1 && iside == 1) {stop("for detecting deterioration (side=\"low\"), QS must a positive numeric value >= 1")}
  else if (QS > 1 && iside == 2) {stop("for detecting improvement, QS must a positive numeric value <= 1")}
  itype <- switch(type, cond = 1, cycl = 2)
  if (is.null(itype)) {
    warning("no valid input, using type=cond (conditional steady-state) as default")
    itype <- 1
  }
  m <- as.integer(m)
  if (is.na(m) || m < 0) {stop("m must a positive integer")}
  .eocusum_adoc_sim(r, k, h, df, coeff, coeff2, QS, iside, itype, m)
}

#' @name eocusum.arloc.h.sim
#' @title Compute alarm threshold of Out of Control EO-CUSUM control charts using simulation
#' @description Compute alarm threshold (Out of Control ARL) of EO-CUSUM control charts using
#'  simulation.

#' @param L0 double. In-Control ARL.
#' @param df DataFrame. First column are Parsonnet Score values within a range of zero to 100 representing
#' the preoperative patient risk. The second column are binary (0/1) outcome values of each operation.
#' @param k A double
#' @param m integer. Number of simulation runs.
#' @param QS double. Defines the performance of a surgeon with the odds ratio ratio of death Q.
#' @param side int 1="upper arm" or 2="lower arm".
#' @param coeff double.
#' @param coeff2 double.
#' @param nc integer. Number of cores.
#' @param OUTPUT boolean. If TRUE verbose output is included, if FALSE a quiet calculation of h is done.
#'
#' @return Returns a single value which is the distance d of the V-mask for a given ARL
#'  and \eqn{\theta}{theta}.
#'
#' @details The function \code{eocusum.arloc.h.sim} determines the control limit for given in-control ARL (L0) by applying a
#' multi-stage search procedure which includes secant rule and the parallel version of \code{\link{eocusum.arloc.sim}}
#' using \code{\link{mclapply}}.
#'
#' @template eocusum.arloc.h.sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum.arloc.h.sim <- function(L0, df, k, coeff, coeff2, m = 100, QS = 1, side = 1, nc = 1, OUTPUT = TRUE) {
  h2 <- 1
  L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arloc.sim, h = h2, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
  if ( OUTPUT ) cat(paste("(i)\t", h2, "\t", L2, "\n"))
  LL <- NULL
  while ( L2 < L0 & h2 < 6 ) {
    L1 <- L2
    h2 <- h2 + 1
    L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arloc.sim, h = h2, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(ii)\t", h2, "\t", L2, "\n"))
    LL <- c(LL, L2)
  }
  if ( L2 < L0 ) {
    x  <- 1:5
    LM <- stats::lm(LL ~ I(x) + I(x^2) )
    beta <- stats::coef(LM)
    p <- beta[2] / beta[3]
    q <- (beta[1] - L0)/beta[3]
    h2 <- -p / 2 + 1*sqrt(p^2 / 4 - q)
    L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arloc.sim, h = h2, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(iii)\t", h2, "\t", L2, "\n"))
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        h2 <- h2 + 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arloc.sim, h = h2, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(iv)a\t", h2, "\t", L2, "\n"))
        }
      h1 <- h2 - 1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        h2 <- h2 - 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arloc.sim, h = h2, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(iv)b\t", h2, "\t", L2, "\n"))
      }
      h1 <- h2 + 1
    }
    } else {
    h1 <- h2 - 1
  }
  h.error <- 1; a.error <- 1; scaling <- 10^3;
  while ( a.error > 1e-4 & h.error > 1e-6 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arloc.sim, h = h3, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(v)\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    a.error <- min( c(abs(L2 - L0), abs(L2 - L1) ) )
    h.error <- abs(h2 - h1)
    if ( h.error < 0.5 / scaling ) {
      if ( L3 < L0 ) {
        h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
        L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arloc.sim, h = h3, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(vi)\t", h3, "\t", L3, "\n"))
      }
      break
    }
  }
  if ( L3 < L0 ) {
    h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
    L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arloc.sim, h = h3, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(vii)\t", h3, "\t", L3, "\n"))
  }
  h <- h3
  h
}

#' @name eocusum.arl.h.sim
#' @title Compute alarm threshold of EO-CUSUM control charts using simulation
#' @description Compute alarm threshold of EO-CUSUM control charts using simulation.
#'
#' @param L0 double. In-Control ARL.
#' @param k A double
#' @param m integer. Number of simulation runs.
#' @param side int 1="upper arm" or 2="lower arm".
#' @param df DataFrame. First column are Parsonnet Score values within a range of zero to 100 representing
#' the preoperative patient risk. The second column are binary (0/1) outcome values of each operation.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta} from the binary
#' logistic regression model. For more information see details.
#' @param yemp boolean. Use emirical outcome value.
#' @param nc integer. Number of cores.
#' @param OUTPUT boolean. If TRUE verbose output is included, if FALSE a quiet calculation of h is done.
#'
#' @return Returns a single value which is the distance d of the V-mask for a given ARL and \eqn{\theta}{theta}.
#'
#' @details The function \code{eocusum.arl.h.sim} determines the control limit for given in-control ARL (L0) by applying a
#' multi-stage search procedure which includes secant rule and the parallel version of \code{\link{eocusum.arl.sim}}
#' using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#' @export
eocusum.arl.h.sim <- function(L0, df, k, coeff, m = 100, yemp = TRUE, side = 1, nc = 1, OUTPUT = TRUE) {
  h2 <- 1
  L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arl.sim, k = k, h = h2, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
  if ( OUTPUT ) cat(paste("(i)\t", h2, "\t", L2, "\n"))
  LL <- NULL
  while ( L2 < L0 & h2 < 6 ) {
    L1 <- L2
    h2 <- h2 + 1
    L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arl.sim, k = k, h = h2, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(ii)\t", h2, "\t", L2, "\n"))
    LL <- c(LL, L2)
  }
  if ( L2 < L0 ) {
    x  <- 1:5
    LM <- stats::lm(LL ~ I(x) + I(x^2) )
    beta <- stats::coef(LM)
    p <- beta[2] / beta[3]
    q <- (beta[1] - L0) / beta[3]
    h2 <- -p / 2 + 1 * sqrt(p^2 / 4 - q)
    L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arl.sim, k = k, h = h2, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(iii)\t", h2, "\t", L2, "\n"))
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        h2 <- h2 + 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arl.sim, k = k, h = h2, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(iv)a\t", h2, "\t", L2, "\n"))
        }
      h1 <- h2 - 1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        h2 <- h2 - 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arl.sim, k = k, h = h2, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(iv)b\t", h2, "\t", L2, "\n"))
      }
      h1 <- h2 + 1
    }
    } else {
    h1 <- h2 - 1
  }
  h.error <- 1; a.error <- 1; scaling <- 10^3;
  while ( a.error > 1e-4 & h.error > 1e-6 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arl.sim, k = k, h = h3, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(v)\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    a.error <- min( c(abs(L2 - L0), abs(L2 - L1) ) )
    h.error <- abs(h2 - h1)
    if ( h.error < 0.5 / scaling ) {
      if ( L3 < L0 ) {
        h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
        L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arl.sim, k = k, h = h3, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(vi)\t", h3, "\t", L3, "\n"))
      }
      break
    }
  }
  if ( L3 < L0 ) {
    h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
    L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum.arl.sim, k = k, h = h3, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(vii)\t", h3, "\t", L3, "\n"))
  }
  h <- h3
  h
}
