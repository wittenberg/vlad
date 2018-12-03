#' @name llr_score
#' @title Compute the log-likelihood ratio score
#' @description Compute the log-likelihood ratio score.
#'
#' @param R0 Double. Odds ratio of death under the null hypotheses.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (\code{0/1})
#'  outcome values of each operation.
#' @param coeff Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model.
#' @param yemp Logical. If \code{TRUE} use observed outcome value, if \code{FALSE} use estimated
#' binary logistc regression model.
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

#' @name bcusum_arl_sim
#' @title Compute ARLs of the Bernoulli CUSUM control charts using simulation
#' @description Compute ARLs of the Bernoulli CUSUM control charts using simulation.
#'
#' @param r Integer Vector. Number of runs.
#' @param h Double. Control Chart limit for detecting deterioration/improvement.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (\code{0/1})
#'  outcome values of each operation.
#' @param R0 Double. Odds ratio of death under the null hypotheses.
#' @param RA Double. Odds ratio of death under the alternative hypotheses.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @author Philipp Wittenberg
#' @export
bcusum_arl_sim <- function(r, h, df, R0 = 1, RA = 2) {
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
  if (is.na(R0) || R0 <= 0) {stop("Odds ratio of death under the null hypotheses 'R0' must a positive numeric value")}
  RA <- as.numeric(RA)
  if (is.na(RA) || RA <= 0) {stop("Odds ratio of death under the alternative hypotheses 'RA' must a positive numeric value")}
  .bcusum_arl_sim(r, h, df, R0, RA)
}

#' @name cusum_arl_sim
#' @title Compute ARLs of Bernoulli CUSUM control charts using simulation
#' @description Compute ARLs of Bernoulli cumulative sum control charts using simulation.
#'
#' @param r Integer Vector. Number of runs.
#' @param h Double. Control Chart limit for detecting deterioration/improvement.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (\code{0/1})
#'  outcome values of each operation.
#' @param R0 Double. Odds ratio of death under the null hypotheses.
#' @param RA Double. Odds ratio of death under the alternative hypotheses.
#'
#' @author Philipp Wittenberg
#'
#' @export
#' @examples
#'\donttest{
#'
#'# This function is deprecated. See bcusum_arl_sim() instead.
#'
#'  }
cusum_arl_sim <- function(r, h, df, R0 = 1, RA = 2) {

  .Deprecated("cusum_arl_sim")
  bcusum_arl_sim(r = r, h = h, df = df, R0 = R0, RA = RA)
}

#' @name racusum_arl_sim
#' @title Compute ARLs of RA-CUSUM control charts using simulation
#' @description Computes the Average Run Length of a risk-adjusted cumulative sum control chart
#'  using simulation.
#'
#' @param r Integer Vector. Number of runs.
#' @inheritParams llr_score
#' @param h Double. Control Chart limit for detecting deterioration/improvement.
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
#' @description Compute Out of Control ARLs of risk-adjusted cumulative sum control charts using
#'  simulation.
#'
#' @inheritParams racusum_arl_sim
#' @param coeff2 Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model of a resampled dataset.
#' @param RQ Double. Defines the performance of a surgeon with the odds ratio ratio of death \code{Q}.
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
  if (is.na(RQ) || RQ <= 0) {stop("RQ must a positive numeric value")}
  .racusum_arloc_sim(r, coeff, coeff2, h, df, R0, RA, RQ)
}

#' @name racusum_adoc_sim
#' @title Compute steady-state ARLs of RA-CUSUM control charts using
#' simulation
#' @description Compute steady-state ARLs of risk-adjusted cumulative sum control charts using
#'  simulation.
#'
#' @inheritParams racusum_arloc_sim
#' @param m Integer. Simulated in-control observations.
#' @param type Character. Default argument is \code{"cond"} for computation of conditional
#' steady-state. Other option is the cyclical steady-state \code{"cycl"}.
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
  if (is.na(RQ) || RQ <= 0) {stop("RQ must a positive numeric value")}
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
#' @description Compute alarm threshold of risk-adjusted cumulative sum control charts using
#'  simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @param R0 Double. Odds ratio of death under the null hypotheses.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#'  in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}.
#'  Detecting improvement in performance with decreased mortality risk by halving the odds ratio of
#'   death \code{RA = 1/2}.
#' @param m Integer. Number of simulation runs.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (\code{0/1})
#'  outcome values of each operation.
#' @param coeff Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model.
#' @param yemp Logical. If \code{TRUE}, use emirical outcome values, else use model.
#' @param nc Integer. Number of cores used for parallel processing.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given in-control ARL.
#'
#' @details The function \code{racusum_arl_h_sim} determines the control limit \code{h} for given
#'  in-control ARL (\code{L0}) by applying a multi-stage search procedure which includes secant
#'  rule and the parallel version of \code{\link{racusum_arl_sim}} using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#'
#' @export
#' @examples
#'\donttest{
#'
#'# This function is deprecated. See racusum_crit_sim() instead.
#'
#'  }
racusum_arl_h_sim <- function(L0, df, coeff, R0 = 1, RA = 2, m = 100, yemp = TRUE, nc = 1, jmax = 4, verbose = FALSE) {

  .Deprecated("racusum_arl_h_sim")
  racusum_crit_sim(L0 = L0, df = df, coeff = coeff, R0 = R0, RA = RA, m = m, yemp = yemp, nc = nc, jmax = jmax, verbose = verbose)
}

#' @name racusum_crit_sim
#' @title Compute alarm threshold of RA-CUSUM control charts using simulation
#' @description Compute alarm threshold of risk-adjusted cumulative sum control charts using
#'  simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @param R0 Double. Odds ratio of death under the null hypotheses.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#'  in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}.
#'  Detecting improvement in performance with decreased mortality risk by halving the odds ratio of
#'   death \code{RA = 1/2}.
#' @param m Integer. Number of simulation runs.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (\code{0/1})
#'  outcome values of each operation.
#' @param coeff Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model.
#' @param yemp Logical. If \code{TRUE}, use emirical outcome values, else use model.
#' @param nc Integer. Number of cores used for parallel processing.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given in-control ARL.
#'
#' @details The function \code{racusum_crit_sim} determines the control limit \code{h} for given
#'  in-control ARL (\code{L0}) by applying a multi-stage search procedure which includes secant
#'  rule and the parallel version of \code{\link{racusum_arl_sim}} using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#'
#' @template racusum_crit_sim
#'
#' @export
racusum_crit_sim <- function(L0, df, coeff, R0 = 1, RA = 2, m = 100, yemp = TRUE, nc = 1, jmax = 4, verbose = FALSE) {
  L0 <- as.integer(L0)
  if (is.na(L0) || L0 <= 0) {stop("Given in-control ARL 'L0' must be a positive integer")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  R0 <- as.numeric(R0)
  if (is.na(R0) || R0 <= 0) {stop("Odds ratio of death under the null hypotheses 'R0' must a positive numeric value")}
  RA <- as.numeric(RA)
  if (is.na(RA) || RA <= 0) {stop("Odds ratio of death under the alternative hypotheses 'RA' must a positive numeric value")}

  for ( h in 1:10 ) {
    L1 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arl_sim, h = h, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      L1 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arl_sim, h = h, df = df, coeff = coeff, R0 = R0, RA = RA, yemp = yemp, mc.cores = nc)))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax

  h
}

#' @name cusum_arl_h_sim
#' @title Compute alarm threshold of Bernoulli CUSUM control charts using simulation
#' @description Compute alarm threshold of Bernoulli cumulative sum control charts using
#' simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @param R0 Double. Odds ratio of death under the null hypotheses.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#'  in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}.
#'  Detecting improvement in performance with decreased mortality risk by halving the odds ratio of
#'  death \code{RA = 1/2}.
#' @param m Integer. Number of simulation runs.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (\code{0/1})
#'  outcome values of each operation.
#' @param nc Integer. Number of cores.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given in-control ARL.
#'
#' @details The function \code{cusum_arl_h_sim} determines the control limit for given in-control
#'  ARL (\code{L0}) by applying a multi-stage search procedure which includes secant rule and the
#'   parallel version of \code{\link{cusum_arl_sim}} using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#'
#' @export
#' @examples
#'\donttest{
#'
#'# This function is deprecated. See bcusum_crit_sim() instead.
#'
#'  }
cusum_arl_h_sim <- function(L0, df, R0 = 1, RA = 2, m = 100, nc = 1, jmax = 4, verbose = FALSE) {

  .Deprecated("cusum_arl_h_sim")
  bcusum_crit_sim(L0 = L0, df = df, R0 = R0, RA = RA, m = m, nc = nc, jmax = jmax, verbose = verbose)
}



#' @name bcusum_crit_sim
#' @title Compute alarm threshold of Bernoulli CUSUM control charts using simulation
#' @description Compute alarm threshold of Bernoulli cumulative sum control charts using
#'  simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @param R0 Double. Odds ratio of death under the null hypotheses.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#'  in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}.
#'  Detecting improvement in performance with decreased mortality risk by halving the odds ratio of
#'  death \code{RA = 1/2}.
#' @param m Integer. Number of simulation runs.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (\code{0/1})
#'  outcome values of each operation.
#' @param nc Integer. Number of cores.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given in-control ARL.
#'
#' @details The function \code{bcusum_crit_sim} determines the control limit for given in-control
#'  ARL (\code{L0}) by applying a multi-stage search procedure which includes secant rule and the
#'   parallel version of \code{\link{bcusum_arl_sim}} using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#' @export
bcusum_crit_sim <- function(L0, df, R0 = 1, RA = 2, m = 100, nc = 1, jmax = 4, verbose = FALSE) {
  L0 <- as.integer(L0)
  if (is.na(L0) || L0 <= 0) {stop("Given in-control ARL 'L0' must be a positive integer")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  R0 <- as.numeric(R0)
  if (is.na(R0) || R0 <= 0) {stop("Odds ratio of death under the null hypotheses 'R0' must a positive numeric value")}
  RA <- as.numeric(RA)
  if (is.na(RA) || RA <= 0) {stop("Odds ratio of death under the alternative hypotheses 'RA' must a positive numeric value")}
  for ( h in 1:10 ) {
    L1 <- mean(do.call(c, parallel::mclapply(1:m, bcusum_arl_sim, h = h, df = df, R0 = R0, RA = RA, mc.cores = nc)))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      L1 <- mean(do.call(c, parallel::mclapply(1:m, bcusum_arl_sim, h = h, df = df, R0 = R0, RA = RA, mc.cores = nc)))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax

  h
}

#' @name racusum_arloc_h_sim
#' @title Compute alarm threshold (Out of Control ARL) of RA-CUSUM control charts using simulation
#' @description Compute alarm threshold (Out of Control ARL) of risk-adjusted cumulative sum
#'  control charts using simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @param R0 Double. Odds ratio of death under the null hypotheses.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#'  in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}.
#'  Detecting improvement in performance with decreased mortality risk by halving the odds ratio of
#'   death \code{RA = 1/2}.
#' @param m Integer. Number of simulation runs.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (\code{0/1})
#'  outcome values of each operation.
#' @param coeff Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model.
#' @param coeff2 Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model of a resampled dataset.
#' @param RQ Double. Defines the performance of a surgeon with the odds ratio ratio of death \code{Q}.
#' @param nc Integer. Number of cores.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given in-control ARL.
#'
#' @details The function \code{racusum_arloc_h_sim} determines the control limit \code{h} for given
#'  in-control ARL (\code{L0}) by applying a multi-stage search procedure which includes secant
#'  rule and the parallel version of \code{\link{racusum_arloc_sim}} using \code{\link{mclapply}}.
#'
#' @template racusum_arloc_h_sim
#'
#' @author Philipp Wittenberg
#' @export
racusum_arloc_h_sim <- function(L0, df, coeff, coeff2, R0 = 1, RA = 2, RQ = 1, m = 100, nc = 1, jmax = 4, verbose = FALSE) {
  L0 <- as.integer(L0)
  if (is.na(L0) || L0 <= 0) {stop("Given in-control ARL 'L0' must be a positive integer")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  R0 <- as.numeric(R0)
  if (is.na(R0) || R0 <= 0) {stop("Odds ratio of death under the null hypotheses 'R0' must a positive numeric value")}
  RA <- as.numeric(RA)
  if (is.na(RA) || RA <= 0) {stop("Odds ratio of death under the alternative hypotheses 'RA' must a positive numeric value")}
  RQ <- as.numeric(RQ)
  if (is.na(RQ) || RQ <= 0) {stop("RQ must a positive numeric value")}
  for ( h in 1:10 ) {
    L1 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      L1 <- mean(do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h, df = df, coeff = coeff, coeff2 = coeff2, R0 = R0, RA = RA, RQ = RQ, mc.cores = nc)))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax

  h
}
