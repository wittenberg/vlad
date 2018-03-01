#' @name llr_score
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
#'  from the binary logistic regression model.
#' @param yemp boolean. If TRUE use observed outcome value, if FALSE use estimated binary logistc
#'  regression model.
#'
#' @return Returns a single value which is the log-likelihood ratio score.
#'
#' @template llr_score
#'
#' @author Philipp Wittenberg
#' @export
llr_score <- function(df, coeff, R0 = 1, RA = 2, yemp = TRUE) {
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff) || length(coeff) != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.na(R0) || R0 <= 0) {stop("Odds ratio of death under the null hypotheses 'R0' must a positive numeric value")}
  R0 <- as.numeric(R0)
  if (is.na(RA) || RA <= 0) {stop("Odds ratio of death under the alternative hypotheses 'RA' must a positive numeric value")}
  RA <- as.numeric(RA)
  if (is.na(yemp) || is.logical(yemp) != "TRUE") {warning("Argument 'yemp' must be logical using TRUE as default value")}
  yemp <- as.logical(yemp)
  .llr_score(df, coeff, R0, RA, yemp)
}

#' @name cusum_arl_sim
#' @title Compute ARLs of the Bernoulli CUSUM control charts using simulation
#' @description Compute ARLs of the Bernoulli CUSUM control charts using simulation.
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
#' @author Philipp Wittenberg
#' @export
cusum_arl_sim <- function(r, h, df, R0 = 1, RA = 2) {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  R0 <- as.numeric(R0)
  if (is.na(R0) || R0 <= 0) {stop("Odds ratio of death under the null hypotheses \"R0\" must a positive numeric value")}
  RA <- as.numeric(RA)
  if (is.na(RA) || RA <= 0) {stop("Odds ratio of death under the alternative hypotheses \"RA\" must a positive numeric value")}
  .cusum_arl_sim(r, h, df, R0, RA)
}

#' @name racusum_arl_sim
#' @title Compute ARLs of RA-CUSUM control charts using simulation
#' @description Compute ARLs of RA-CUSUM control charts using simulation.
#'
#' @param r Integer vector. Number of runs.
#' @inheritParams llr_score
#' @param h double. Control Chart limit for detecting deterioration/improvement.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template racusum_arl_sim
#'
#' @author Philipp Wittenberg
#' @export
racusum_arl_sim <- function(r, coeff, h, df, R0 = 1, RA = 2, yemp = TRUE) {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  if (is.null(coeff) || is.na(coeff) || length(coeff) != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  R0 <- as.numeric(R0)
  if (is.na(R0) || R0 <= 0) {stop("Odds ratio of death under the null hypotheses 'R0' must a positive numeric value")}
  RA <- as.numeric(RA)
  if (is.na(RA) || RA <= 0) {stop("Odds ratio of death under the alternative hypotheses 'RA' must a positive numeric value")}
  if (is.na(yemp) || is.logical(yemp) != "TRUE") {warning("Argument 'yemp' must be logical using TRUE as default value")}
  yemp <- as.logical(yemp)
  .racusum_arl_sim(r, coeff, h, df, R0, RA, yemp)
}

#' @name racusum_arloc_sim
#' @title Compute Out of Control ARLs of RA-CUSUM control charts using simulation
#' @description Compute Out of Control ARLs of RA-CUSUM control charts using simulation.
#'
#' @inheritParams racusum_arl_sim
#' @param coeff2 NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model of a resampled dataset.
#' @param RQ double. Defines the performance of a surgeon with the odds ratio ratio of death Q.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template racusum_arloc_sim
#'
#' @author Philipp Wittenberg
#'
#' @export
racusum_arloc_sim <- function(r, coeff, coeff2, h, df, R0 = 1, RA = 2, RQ = 1) {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  if (is.null(coeff) || is.na(coeff) || length(coeff) != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.null(coeff2) || is.na(coeff2) || length(coeff2) != 2) {stop("Model coefficients 'coeff2' must be a numeric vector with two elements")}
  coeff2 <- as.vector(coeff2)
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.na(R0) || R0 <= 0) {stop("Odds ratio of death under the null hypotheses 'R0' must a positive numeric value")}
  RA <- as.numeric(RA)
  if (is.na(RA) || RA <= 0) {stop("Odds ratio of death under the alternative hypotheses 'RA' must a positive numeric value")}
  RQ <- as.numeric(RQ)
  if (is.na(RQ) || RQ < 0) {stop("RQ must a positive numeric value")}
  .racusum_arloc_sim(r, coeff, coeff2, h, df, R0, RA, RQ)
}

#' @name racusum_adoc_sim
#' @title Compute steady-state ARLs of RA-CUSUM control charts using
#' simulation
#' @description Compute steady-state ARLs of RA-CUSUM control charts using simulation.
#'
#' @inheritParams racusum_arloc_sim
#' @param m Integer. Simulated in-control observations.
#' @param type character. Default argument is "cond" for computation of conditional steady-state.
#' Other option is the cyclical steady-state "cycl".
#' @return Returns a single value which is the Run Length.
#'
#' @template racusum_adoc_sim
#'
#' @author Philipp Wittenberg
#'
#' @export
racusum_adoc_sim <- function(r, coeff, coeff2, h, df, R0 = 1, RA = 2, RQ = 1, m = 50, type = "cond") {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  if (is.null(coeff) || is.na(coeff) || length(coeff) != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.null(coeff2) || is.na(coeff2) || length(coeff2) != 2) {stop("Model coefficients 'coeff2' must be a numeric vector with two elements")}
  coeff2 <- as.vector(coeff2)
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  R0 <- as.numeric(R0)
  if (is.na(R0) || R0 <= 0) {stop("Odds ratio of death under the null hypotheses 'R0' must a positive numeric value")}
  RA <- as.numeric(RA)
  if (is.na(RA) || RA <= 0) {stop("Odds ratio of death under the alternative hypotheses 'RA' must a positive numeric value")}
  RQ <- as.numeric(RQ)
  if (is.na(RQ) || RQ < 0) {stop("RQ must a positive numeric value")}
  m <- as.integer(m)
  if (is.na(m) || m < 0) {stop("m must a positive integer")}
  itype <- switch(type, cond = 1, cycl = 2)
  if (is.null(itype)) {
    warning("no valid input, using type=cond (conditional steady-state) as default")
    itype <- 1
  }
  .racusum_adoc_sim(r, coeff, coeff2, h, df, R0, RA, RQ, m, itype)
}

#' @name racusum_arl_h_sim
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
#'  from the binary logistic regression model.
#' @param yemp boolean. If TRUE, use emirical outcome values, else use model.
#' @param nc integer. Number of cores used for parallel processing.
#' @param verbose boolean. If TRUE verbose output is included, if FALSE a quiet calculation of h is done.
#'
#' @return Returns a single value which is the control limit h for a given in-control ARL.
#'
#' @details The function \code{racusum_arl_h_sim} determines the control limit for given in-control ARL (L0) by applying a
#' multi-stage search procedure which includes secant rule and the parallel version of \code{\link{racusum_arl_sim}}
#' using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#'
#' @template racusum_arl_h_sim
#'
#' @export
racusum_arl_h_sim <- function(L0, df, coeff, R0 = 1, RA = 2, m = 100, yemp = TRUE, nc = 1, verbose = FALSE) {
  h2 <- 1
  L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arl_sim, h = h2, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
  if ( verbose ) cat(paste("(i)\t", h2, "\t", L2, "\n"))
  LL <- NULL
  while ( L2 < L0 & h2 < 6 ) {
    L1 <- L2
    h2 <- h2 + 1
    L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arl_sim, h = h2, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
    if ( verbose ) cat(paste("(ii)\t", h2, "\t", L2, "\n"))
    LL <- c(LL, L2)
  }
  if ( L2 < L0 ) {
    x  <- 1:5
    LM <- stats::lm(LL ~ I(x) + I(x^2) )
    beta <- stats::coef(LM)
    p <- beta[2] / beta[3]
    q <- (beta[1] - L0) / beta[3]
    h2 <- -p / 2 + 1 * sqrt(p^2 / 4 - q)
    L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arl_sim, h = h2, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
    if ( verbose ) cat(paste("(iii)\t", h2, "\t", L2, "\n"))
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        h2 <- h2 + 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arl_sim, h = h2, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
        if ( verbose ) cat(paste("(iv)a\t", h2, "\t", L2, "\n"))
        }
      h1 <- h2 - 1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        h2 <- h2 - 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arl_sim, h = h2, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
        if ( verbose ) cat(paste("(iv)b\t", h2, "\t", L2, "\n"))
      }
      h1 <- h2 + 1
    }
    } else {
    h1 <- h2 - 1
  }
  h.error <- 1
  a.error <- 1
  scaling <- 10^3
  while ( a.error > 1e-4 & h.error > 1e-6 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arl_sim, h = h3, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
    if ( verbose ) cat(paste("(v)\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    a.error <- min( c(abs(L2 - L0), abs(L2 - L1) ) )
    h.error <- abs(h2 - h1)
    if ( h.error < 0.5 / scaling ) {
      if ( L3 < L0 ) {
        h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
        L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arl_sim, h = h3, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
        if ( verbose ) cat(paste("(vi)\t", h3, "\t", L3, "\n"))
      }
      break
    }
  }
  if ( L3 < L0 ) {
    h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
    L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arl_sim, h = h3, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
    if ( verbose ) cat(paste("(vii)\t", h3, "\t", L3, "\n"))
  }
  h <- h3
  h
}

#' @name cusum_arl_h_sim
#' @title Compute alarm threshold of the Bernoulli CUSUM control charts using simulation
#' @description Compute alarm threshold of the Bernoulli CUSUM control charts using simulation.
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
#' @param nc integer. Number of cores.
#' @param verbose boolean. If TRUE verbose output is included, if FALSE a quiet calculation of h is done.
#'
#' @return Returns a single value which is the control limit h for a given in-control ARL.
#'
#' @details The function \code{cusum_arl_h_sim} determines the control limit for given in-control ARL (L0) by applying a
#' multi-stage search procedure which includes secant rule and the parallel version of \code{\link{cusum_arl_sim}}
#' using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#' @export
cusum_arl_h_sim <- function(L0, df, R0 = 1, RA = 2, m = 100, nc = 1, verbose = FALSE) {
  h2 <- 1
  L2 <- mean(do.call(c, parallel::mclapply(1:m, cusum_arl_sim, h = h2, df = df, R0 = R0, RA = RA, mc.cores = nc)))
  if ( verbose ) cat(paste("(i)\t", h2, "\t", L2, "\n"))
  LL <- NULL
  while ( L2 < L0 & h2 < 6 ) {
    L1 <- L2
    h2 <- h2 + 1
    L2 <- mean(do.call(c, parallel::mclapply(1:m, cusum_arl_sim, h = h2, df = df, R0 = R0, RA = RA, mc.cores = nc)))
    if ( verbose ) cat(paste("(ii)\t", h2, "\t", L2, "\n"))
    LL <- c(LL, L2)
  }
  if ( L2 < L0 ) {
    x  <- 1:5
    LM <- stats::lm(LL ~ I(x) + I(x^2) )
    beta <- stats::coef(LM)
    p <- beta[2] / beta[3]
    q <- (beta[1] - L0) / beta[3]
    h2 <- -p / 2 + 1 * sqrt(p^2 / 4 - q)
    L2 <- mean(do.call(c, parallel::mclapply(1:m, cusum_arl_sim, h = h2, df = df, R0 = R0, RA = RA, mc.cores = nc)))
    if ( verbose ) cat(paste("(iii)\t", h2, "\t", L2, "\n"))
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        h2 <- h2 + 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, cusum_arl_sim, h = h2, df = df, R0 = R0, RA = RA, mc.cores = nc)))
        if ( verbose ) cat(paste("(iv)a\t", h2, "\t", L2, "\n"))
        }
      h1 <- h2 - 1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        h2 <- h2 - 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, cusum_arl_sim, h = h2, df = df, R0 = R0, RA = RA, mc.cores = nc)))
        if ( verbose ) cat(paste("(iv)b\t", h2, "\t", L2, "\n"))
      }
      h1 <- h2 + 1
    }
  } else {
    h1 <- h2 - 1
  }
  h.error <- 1
  a.error <- 1
  scaling <- 10^3
  while ( a.error > 1e-4 & h.error > 1e-6 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- mean(do.call(c, parallel::mclapply(1:m, cusum_arl_sim, h = h3, df = df, R0 = R0, RA = RA, mc.cores = nc)))
    if ( verbose ) cat(paste("(v)\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    a.error <- min( c(abs(L2 - L0), abs(L2 - L1) ) )
    h.error <- abs(h2 - h1)
    if ( h.error < 0.5 / scaling ) {
      if ( L3 < L0 ) {
        h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
        L3 <- mean(do.call(c, parallel::mclapply(1:m, cusum_arl_sim, h = h3, df = df, R0 = R0, RA = RA, mc.cores = nc)))
        if ( verbose ) cat(paste("(vi)\t", h3, "\t", L3, "\n"))
      }
      break
    }
  }
  if ( L3 < L0 ) {
    h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
    L3 <- mean(do.call(c, parallel::mclapply(1:m, cusum_arl_sim, h = h3, df = df, R0 = R0, RA = RA, mc.cores = nc)))
    if ( verbose ) cat(paste("(vii)\t", h3, "\t", L3, "\n"))
  }
  h <- h3
  h
}

#' @name racusum_arloc_h_sim
#' @title Compute alarm threshold (Out of Control ARL) of RA-CUSUM control charts using simulation
#' @description Compute alarm threshold (Out of Control ARL) of RA-CUSUM control charts using
#' simulation.
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
#'  from the binary logistic regression model.
#' @param coeff2 NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model of a resampled dataset.
#' @param RQ double. Defines the performance of a surgeon with the odds ratio ratio of death Q.
#' @param nc integer. Number of cores.
#' @param verbose boolean. If TRUE verbose output is included, if FALSE a quiet calculation of h is done.
#'
#' @return Returns a single value which is the control limit h for a given in-control ARL.
#'
#' @details The function \code{racusum_arloc_h_sim} determines the control limit for given in-control ARL (L0) by applying a
#' multi-stage search procedure which includes secant rule and the parallel version of \code{\link{racusum_arloc_sim}}
#' using \code{\link{mclapply}}.
#'
#' @template racusum_arloc_h_sim
#'
#' @author Philipp Wittenberg
#' @export
racusum_arloc_h_sim <- function(L0, df, coeff, coeff2, R0 = 1, RA = 2, RQ = 1, m = 100, nc = 1, verbose = FALSE) {
  h2 <- 1
  L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h2, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
  if ( verbose ) cat(paste("(i)\t", h2, "\t", L2, "\n"))
  LL <- NULL
  while ( L2 < L0 & h2 < 6 ) {
    L1 <- L2
    h2 <- h2 + 1
    L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h2, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
    if ( verbose ) cat(paste("(ii)\t", h2, "\t", L2, "\n"))
    LL <- c(LL, L2)
  }
  if ( L2 < L0 ) {
    x  <- 1:5
    LM <- stats::lm(LL ~ I(x) + I(x^2) )
    beta <- stats::coef(LM)
    p <- beta[2] / beta[3]
    q <- (beta[1] - L0) / beta[3]
    h2 <- -p / 2 + 1 * sqrt(p^2 / 4 - q)
    L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h2, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
    if ( verbose ) cat(paste("(iii)\t", h2, "\t", L2, "\n"))
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        h2 <- h2 + 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h2, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
        if ( verbose ) cat(paste("(iv)a\t", h2, "\t", L2, "\n"))
        }
      h1 <- h2 - 1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        h2 <- h2 - 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h2, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
        if ( verbose ) cat(paste("(iv)b\t", h2, "\t", L2, "\n"))
      }
      h1 <- h2 + 1
    }
    } else {
    h1 <- h2 - 1
  }
  h.error <- 1
  a.error <- 1
  scaling <- 10^3
  while ( a.error > 1e-4 & h.error > 1e-6 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h3, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
    if ( verbose ) cat(paste("(v)\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    a.error <- min( c(abs(L2 - L0), abs(L2 - L1) ) )
    h.error <- abs(h2 - h1)
    if ( h.error < 0.5 / scaling ) {
      if ( L3 < L0 ) {
        h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
        L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h3, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
        if ( verbose ) cat(paste("(vi)\t", h3, "\t", L3, "\n"))
      }
      break
    }
  }
  if ( L3 < L0 ) {
    h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
    L3 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h3, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
    if ( verbose ) cat(paste("(vii)\t", h3, "\t", L3, "\n"))
  }
  h <- h3
  h
}
