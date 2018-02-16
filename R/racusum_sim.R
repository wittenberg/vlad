#' @name loglikelihood
#' @title Compute the log-likelihood ratio score
#' @description Compute the log-likelihood ratio score.
#'
#' @param R0 double. Odds ratio of death under the null hypotheses.
#' @param RA double. Odds ratio of death under the alternative hypotheses.
#'  Detecting deterioration in performance with increased mortality risk by doubling the odds Ratio
#'  RA=2. Detecting improvement in performance with decreased mortality risk by halving the odds
#'  ratio of death RA=1/2.
#' @param df DataFrame. First column are Parsonnet Score values within a range of zero to 100 representing
#' the preoperative patient risk. The second column are binary (0/1) outcome values of each operation.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model. For more information see details.
#' @param yemp boolean. If TRUE use observed outcome value, if FALSE use estimated binary logistc
#'  regression model.
#'
#' @return Returns a single value which is the log-likelihood ratio score.
#'
#' @details Using .... include equation of log-likelihood ratio here and use
#'  \code{\link{gettherisk}} to calulate the risk of failure.
#'
#' @template loglikelihood
#'
#' @author Philipp Wittenberg
#' @export
loglikelihood <- function(df, coeff, R0 = 1, RA = 2, yemp = TRUE) {
  .loglikelihood(
    as.data.frame(df),
    as.vector(coeff),
    as.numeric(R0),
    as.numeric(RA),
    as.logical(yemp)
  )
}

#' @name racusum.arl.nonRA.sim
#' @title Compute ARLs of Non RA-CUSUM control charts using simulation
#' @description Compute ARLs of Non RA-CUSUM control charts using simulation.
#'
#' @param r Integer vector. Number of runs.
#' @param h double. Control Chart limit for detecting deterioration/improvement.
#' @param df DataFrame. First column are Parsonnet Score values within a range of zero to 100 representing
#' the preoperative patient risk. The second column are binary (0/1) outcome values of each operation.
#' @param R0 double. Odds ratio of death under the null hypotheses.
#' @param RA double. Odds ratio of death under the alternative hypotheses.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @details Describe the algorithm for calulating in control run length ...
#'
#' @references \insertRef{Steiner.etal_2000}{VLAD2}
#'
#' @author Philipp Wittenberg
#' @export
racusum.arl.nonRA.sim <- function(r, h, df, R0 = 1, RA = 2) {
  .racusum_arl_nonRA(
    as.integer(r),
    as.numeric(h),
    as.data.frame(df),
    as.numeric(R0),
    as.numeric(RA)
  )
}

#' @name racusum.arl.sim
#' @title Compute ARLs of RA-CUSUM control charts using simulation
#' @description Compute ARLs of RA-CUSUM control charts using simulation.
#'
#' @param r Integer vector. Number of runs.
#' @inheritParams loglikelihood
#' @param h double. Control Chart limit for detecting deterioration/improvement.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @details Describe the algorithm for calulating in control run length ...
#'
#' @template racusum.arl.sim
#'
#' @author Philipp Wittenberg
#' @export
racusum.arl.sim <- function(r, coeff, h, df, R0 = 1, RA = 2, yemp = TRUE) {
  .racusum_arl_sim(
    as.integer(r),
    as.vector(coeff),
    as.numeric(h),
    as.data.frame(df),
    as.numeric(R0),
    as.numeric(RA),
    as.logical(yemp)
  )
}

#' @name racusum.arloc.sim
#' @title Compute Out of Control ARLs of RA-CUSUM control charts using simulation
#' @description Compute Out of Control ARLs of RA-CUSUM control charts using simulation.
#'
#' @inheritParams racusum.arl.sim
#' @param coeff2 NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model of a resampled dataset. For more information see
#'  details.
#' @param RQ double. Defines the performance of a surgeon with the odds ratio ratio of death Q.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @details Describe the resampling algorithm for calulating out of control run length.
#'
#' @template racusum.arloc.sim
#'
#' @author Philipp Wittenberg
#'
#' @export
racusum.arloc.sim <- function(r, coeff, coeff2, h, df, R0 = 1, RA = 2, RQ = 1) {
  .racusum_arloc_sim(
    as.integer(r),
    as.vector(coeff),
    as.vector(coeff2),
    as.numeric(h),
    as.data.frame(df),
    as.numeric(R0),
    as.numeric(RA),
    as.numeric(RQ)
  )
}

#' @name racusum.adoc.sim
#' @title Compute conditional steady-state ARLs of RA-CUSUM control charts using
#' simulation
#' @description Compute conditional steady-state ARLs of RA-CUSUM control charts using simulation.
#'
#' @inheritParams racusum.arloc.sim
#' @param m Integer. Simulated in-control observations.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @details Describe the resampling algorithm for calulating out of control run length.
#'
#' @template racusum.adoc.sim
#'
#' @author Philipp Wittenberg
#'
#' @export
racusum.adoc.sim <- function(r, coeff, coeff2, h, df, R0 = 1, RA = 2, RQ = 1, m = 50) {
  .racusum_adoc_sim(
    as.integer(r),
    as.vector(coeff),
    as.vector(coeff2),
    as.numeric(h),
    as.data.frame(df),
    as.numeric(R0),
    as.numeric(RA),
    as.numeric(RQ),
    as.integer(m)
  )
}

#' @name racusum.adoc2.sim
#' @title Compute cycliclal steady-state ARLs of RA-CUSUM control charts using simulation
#' @description Compute cycliclal steady-state ARLs of RA-CUSUM control charts using simulation.
#'
#' @inheritParams racusum.arloc.sim
#' @param m Integer. Simulated in-control observations.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @details Describe the resampling algorithm for calulating out of control run length.
#'
#' @template racusum.adoc2.sim
#'
#' @author Philipp Wittenberg
#' @export
racusum.adoc2.sim <- function(r, coeff, coeff2, h, df, R0 = 1, RA = 2, RQ = 1, m = 50) {
  .racusum_adoc2_sim(
    as.integer(r),
    as.vector(coeff),
    as.vector(coeff2),
    as.numeric(h),
    as.data.frame(df),
    as.numeric(R0),
    as.numeric(RA),
    as.numeric(RQ),
    as.integer(m)
  )
}

#' @name racusum.arl.h.sim
#' @title Compute alarm threshold of RA-CUSUM control charts using simulation
#' @description Compute alarm threshold of RA-CUSUM control charts using simulation.
#'
#' @param L0 double. Prespecified in-control Average Run Length.
#' @param R0 double. Odds ratio of death under the null hypotheses.
#' @param RA double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#'  in performance with increased mortality risk by doubling the odds Ratio RA=2. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  RA=1/2.
#' @param m integer. Number of simulation runs.
#' @param df DataFrame. First column are Parsonnet Score values within a range of zero to 100 representing
#' the preoperative patient risk. The second column are binary (0/1) outcome values of each operation.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model. For more information see details.
#' @param yemp boolean. If TRUE, use emirical outcome values, else use model.
#' @param nc integer. Number of cores used for parallel processing.
#' @param OUTPUT boolean. If TRUE verbose output is included, if FALSE a quiet calculation of h is done.
#'
#' @return Returns a single value which is the control limit h for a given in-control ARL.
#'
#' @details The function \code{racusum.arl.h.sim} determines the control limit for given in-control ARL (L0) by applying a
#' multi-stage search procedure which includes secant rule and the parallel version of \code{\link{racusum.arl.sim}}
#' using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#' @examples
#' \dontrun{
#' library("vlad")
#' library("spcadjust")
#' data("cardiacsurgery")
#' cardiacsurgery <- dplyr::mutate(cardiacsurgery, phase=factor(ifelse(date < 2*365, "I", "II")))
#' S2 <- subset(cardiacsurgery, c(surgeon==2), c("phase", "Parsonnet", "status"))
#' # subset phase I (In-control) of surgeons 2
#' S2I <- subset(S2, c(phase=="I"))
#' # estimate coefficients from logit model
#' coeff1 <- round(coef(glm(status~Parsonnet, data=S2I, family="binomial")), 3)
#'
#' racusum.arl.h.sim(L0=740, df=S2I, coeff=coeff1, m=10^2, nc=4)
#'}
#' @export
racusum.arl.h.sim <- function(L0, df, coeff, R0 = 1, RA = 2, m = 100, yemp = TRUE, nc = 1, OUTPUT = TRUE) {
  h2 <- 1
  L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.sim, h = h2, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
  if ( OUTPUT ) cat(paste("(i)\t", h2, "\t", L2, "\n"))
  LL <- NULL
  while ( L2 < L0 & h2 < 6 ) {
    L1 <- L2
    h2 <- h2 + 1
    L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.sim, h = h2, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
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
    L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.sim, h = h2, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(iii)\t", h2, "\t", L2, "\n"))
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        h2 <- h2 + 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.sim, h = h2, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(iv)a\t", h2, "\t", L2, "\n"))
        }
      h1 <- h2 - 1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        h2 <- h2 - 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.sim, h = h2, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
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
    L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.sim, h = h3, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
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
        L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.sim, h = h3, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(vi)\t", h3, "\t", L3, "\n"))
      }
      break
    }
  }
  if ( L3 < L0 ) {
    h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
    L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.sim, h = h3, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(vii)\t", h3, "\t", L3, "\n"))
  }
  h <- h3
  h
}

#' @name racusum.arl.nonRA.h.sim
#' @title Compute alarm threshold of Non RA-CUSUM control charts using simulation
#' @description Compute alarm threshold of Non RA-CUSUM control charts using simulation.
#'
#' @param L0 double. Prespecified Average Run Length.
#' @param R0 double. Odds ratio of death under the null hypotheses.
#' @param RA double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#'  in performance with increased mortality risk by doubling the odds Ratio RA=2. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  RA=1/2.
#' @param m integer. Number of simulation runs.
#' @param df DataFrame. First column are Parsonnet Score values within a range of zero to 100 representing
#' the preoperative patient risk. The second column are binary (0/1) outcome values of each operation.
#' @param nc integer. Number of cores.
#' @param sv double. Starting value of secant method.
#' @param stepsize double. Parameter for tuning  increase/decrease of secant method.
#' @param OUTPUT swrgvares
#'
#' @return Returns a single value which is the control limit h for a given in-control ARL.
#'
#' @details The function \code{racusum.arl.nonRA.h.sim} determines the control limit for given in-control ARL (L0) by applying a
#' multi-stage search procedure which includes secant rule and the parallel version of \code{\link{racusum.arl.nonRA.sim}}
#' using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#' @export
racusum.arl.nonRA.h.sim <- function(L0, df, R0 = 1, RA = 2, m = 100, nc = 1, sv = 0.5, stepsize = 0.5, OUTPUT = TRUE) {
  h2 <- 1
  L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.nonRA.sim, h = h2, df = df, R0 = R0, RA = RA, mc.cores = nc)))
  if ( OUTPUT ) cat(paste("(i)\t", h2, "\t", L2, "\n"))
  LL <- NULL
  while ( L2 < L0 & h2 < 6 ) {
    L1 <- L2
    h2 <- h2 + 1
    L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.nonRA.sim, h = h2, df = df, R0 = R0, RA = RA, mc.cores = nc)))
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
    L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.nonRA.sim, h = h2, df = df, R0 = R0, RA = RA, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(iii)\t", h2, "\t", L2, "\n"))
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        h2 <- h2 + 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.nonRA.sim, h = h2, df = df, R0 = R0, RA = RA, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(iv)a\t", h2, "\t", L2, "\n"))
        }
      h1 <- h2 - 1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        h2 <- h2 - 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.nonRA.sim, h = h2, df = df, R0 = R0, RA = RA, mc.cores = nc)))
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
    L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.nonRA.sim, h = h3, df = df, R0 = R0, RA = RA, mc.cores = nc)))
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
        L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.nonRA.sim, h = h3, df = df, R0 = R0, RA = RA, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(vi)\t", h3, "\t", L3, "\n"))
      }
      break
    }
  }
  if ( L3 < L0 ) {
    h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
    L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arl.nonRA.sim, h = h3, df = df, R0 = R0, RA = RA, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(vii)\t", h3, "\t", L3, "\n"))
  }
  h <- h3
  h
}

#' @name racusum.arloc.h.sim
#' @title Compute alarm threshold (Out of Control ARL) of RA-CUSUM control charts using simulation
#' @description Compute alarm threshold (Out of Control ARL) of RA-CUSUM control charts using
#' simulation.
#'
#' @param L0 double. Prespecified Average Run Length.
#' @param R0 double. Odds ratio of death under the null hypotheses.
#' @param RA double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#'  in performance with increased mortality risk by doubling the odds Ratio RA=2. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  RA=1/2.
#' @param m integer. Number of simulation runs.
#' @param df DataFrame. First column are Parsonnet Score values within a range of zero to 100 representing
#' the preoperative patient risk. The second column are binary (0/1) outcome values of each operation.
#' @param coeff NumericVector
#' @param coeff2 NumericVector
#' @param RQ double. Defines the performance of a surgeon with the odds ratio ratio of death Q.
#' @param nc integer. Number of cores.
#' @param OUTPUT boolean. If TRUE verbose output is included, if FALSE a quiet calculation of h is done.
#'
#' @return Returns a single value which is the control limit h for a given in-control ARL.
#'
#' @details The function \code{racusum.arloc.h.sim} determines the control limit for given in-control ARL (L0) by applying a
#' multi-stage search procedure which includes secant rule and the parallel version of \code{\link{racusum.arloc.sim}}
#' using \code{\link{mclapply}}.
#'
#' @template racusum.arloc.h.sim
#'
#' @author Philipp Wittenberg
#' @export
racusum.arloc.h.sim <- function(L0, df, coeff, coeff2, R0 = 1, RA = 2, RQ = 1, m = 100, nc = 1, OUTPUT = TRUE) {
  h2 <- 1
  L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arloc.sim, h = h2, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
  if ( OUTPUT ) cat(paste("(i)\t", h2, "\t", L2, "\n"))
  LL <- NULL
  while ( L2 < L0 & h2 < 6 ) {
    L1 <- L2
    h2 <- h2 + 1
    L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arloc.sim, h = h2, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
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
    L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arloc.sim, h = h2, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(iii)\t", h2, "\t", L2, "\n"))
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        h2 <- h2 + 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arloc.sim, h = h2, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(iv)a\t", h2, "\t", L2, "\n"))
        }
      h1 <- h2 - 1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        h2 <- h2 - 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arloc.sim, h = h2, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
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
    L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arloc.sim, h = h3, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
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
        L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arloc.sim, h = h3, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
        if ( OUTPUT ) cat(paste("(vi)\t", h3, "\t", L3, "\n"))
      }
      break
    }
  }
  if ( L3 < L0 ) {
    h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
    L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum.arloc.sim, h = h3, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
    if ( OUTPUT ) cat(paste("(vii)\t", h3, "\t", L3, "\n"))
  }
  h <- h3
  h
}
