#' @name optimal_k
#' @title Compute optimal k
#' @description Compute optimal k.
#'
#' @param QA Double. Defines the performance of a surgeon with the odds ratio ratio of death
#' \code{Q}.
#' @param df Data Frame. First column Parsonnet Score and second column outcome of each operation.
#' @param coeff Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model. For more information see details.
#' @param yemp Logical. If \code{TRUE} use observed outcome value, if \code{FALSE} use estimated
#' binary logistc regression model.
#'
#' @return Returns a single value which is the approximate optimal \code{k} for a set of given
#' Parsonnet scores.
#'
#' @details Formula deterioration:  \deqn{ k{det} = \frac{Q{A} - 1 - log(Q{A})}{log(Q{A})}\bar{p} , Q{A} > 1    }
#'          Formula improvement:    \deqn{ k{imp} = \frac{1 - Q{A} + log(Q{A})}{log(Q{A})}\bar{p} , Q{A} < 1    }
#'
#' @template optimal_k
#'
#' @author Philipp Wittenberg
#' @export
optimal_k <- function(QA, df, coeff, yemp = TRUE) {
  if (is.null(QA) || is.na(QA) || QA <= 0) {stop("QA must a positive numeric value")}
  QA <- as.numeric(QA)
  #if (is.null(parsonnetscores) || is.na(parsonnetscores) || is.vector(parsonnetscores) != "TRUE") {stop("Argument 'parsonnetscore' must be an integer value")}
  #parsonnetscores <- as.vector(parsonnetscores)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.na(yemp) || is.logical(yemp) != "TRUE") {warning("Argument 'yemp' must be logical using TRUE as default value")}
  yemp <- as.logical(yemp)
  .optimal_k(QA, df, coeff, yemp)
}

#' @name gettherisk
#' @title Compute Risk of death
#' @description Compute Risk of death.
#'
#' @param parsonnetscore Integer. Parsonnet Score.
#' @param coeff Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model.
#'
#' @return Returns a single value which is the expected risk based on a risk model.
#'
#' @template gettherisk
#'
#' @author Philipp Wittenberg
#'
#' @export
gettherisk <- function(parsonnetscore, coeff) {
  if (is.null(parsonnetscore) || is.na(parsonnetscore) || is.integer(parsonnetscore) != "TRUE") {stop("Argument 'parsonnetscore' must be an integer value")}
  parsonnetscore <- as.integer(parsonnetscore)
  if (is.null(coeff) || is.na(coeff) || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  .gettherisk(parsonnetscore, coeff)
}

#' @name calceo
#' @title Compute Expected minus Observed value
#' @description Compute Expected minus Observed value.
#'
#' @param df Data Frame. First column Parsonnet Score and second column outcome of each operation.
#' @param coeff Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model.
#' @param yemp Logical. If \code{TRUE} use observed outcome value, if \code{FALSE} use estimated
#' binary logistc regression model.
#'
#' @return Returns a single value which is the difference between expected risk and observed
#' outcome.
#'
#' @template calceo
#'
#' @author Philipp Wittenberg
#' @export
calceo <- function(df, coeff, yemp = TRUE) {
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff) || length(coeff) != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.na(yemp) || is.logical(yemp) != "TRUE") {warning("Argument 'yemp' must be logical using TRUE as default value")}
  yemp <- as.logical(yemp)
  .calceo(df, coeff, yemp)
}

#' @name eocusum_arl_sim
#' @title Compute ARLs of EO-CUSUM control charts using simulation
#' @description Compute ARLs of EO-CUSUM control charts using simulation.
#'
#' @param r Integer. Number of of simulation runs.
#' @param k Double. Reference value of the CUSUM control chart. Either \code{0} or a positive
#' value. Can be determined with function \code{\link{optimal_k}}.
#' @param h Double. Decision interval (alarm limit, threshold) of the CUSUM control chart.
#' @inheritParams calceo
#' @param side Character. Default is \code{"low"} to calculate ARL for the upper arm of the V-mask.
#'  If side = \code{"up"}, calculate the lower arm of the V-mask.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template eocusum_arl_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_arl_sim <- function(r, k, h, df, coeff, yemp = TRUE, side = "low") {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff) || length(coeff) != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.na(yemp) || is.logical(yemp) != "TRUE") {warning("Argument 'yemp' must be logical using TRUE as default value")}
  yemp <- as.logical(yemp)
  iside <- switch(as.character(side), low = 1, up = 2)
  if (is.null(iside)) {
    warning("No valid input, using side='low' (deterioration) as default")
    iside <- 1
  }
  .eocusum_arl_sim(r, k, h, df, coeff, yemp, iside)
}

#' @name eocusum_arloc_sim
#' @title Compute Out of Control ARLs of EO-CUSUM control charts using simulation
#' @description Compute Out of Control ARLs of EO-CUSUM control charts using simulation.
#'
#' @inheritParams eocusum_arl_sim
#' @param coeff2 Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model of a resampled dataset.
#' @param QS Double. Defines the performance of a surgeon with the odds ratio ratio of death
#' \code{Q}.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template eocusum_arloc_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_arloc_sim <- function(r, k, h, df, coeff, coeff2, QS = 1, side = "low") {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.null(coeff2) || is.na(coeff2) || length(coeff2) != 2) {stop("Model coefficients 'coeff2' must be a numeric vector with two elements")}
  coeff2 <- as.vector(coeff2)
  iside <- switch(as.character(side), low = 1, up = 2)
  if (is.null(iside)) {
    warning("No valid input, using side='low' (deterioration) as default")
    iside <- 1
  }
  QS <- as.numeric(QS)
  if (is.na(QS) || QS <= 0) {stop("QS must a positive numeric value")}
  else if (QS < 1 && iside == 1) {stop("For detecting deterioration (side='low') QS must a positive numeric value >= 1")}
  else if (QS > 1 && iside == 2) {stop("For detecting improvement (side='up') QS must a positive numeric value <= 1")}
  .eocusum_arloc_sim(r, k, h, df, coeff, coeff2, QS, iside)
}

#' @name eocusum_adoc_sim
#' @title Compute steady-state ARLs of EO-CUSUM control charts using simulation
#' @description Compute steady-state ARLs of EO-CUSUM control charts using simulation.
#'
#' @inheritParams eocusum_arloc_sim
#' @param m Integer. Simulated in-control observations.
#' @param type Character. Default argument is \code{"cond"} for computation of conditional
#' steady-state. Other option is the cyclical steady-state \code{"cycl"}.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @author Philipp Wittenberg
#'
#' @export
#' @examples
#'\donttest{
#'
#'# This function is deprecated. See eocusum_ad_sim() instead.
#'
#'  }
eocusum_adoc_sim <- function(r, k, h, df, coeff, coeff2, QS = 1, side = "low", type = "cond", m = 50) {

  .Deprecated("eocusum_adoc_sim")
  eocusum_ad_sim(r = r, k = k, h = h, df = df, coeff = coeff, coeff2 = coeff2, QS = QS, side = side, type = type, m = m)
}

#' @name eocusum_ad_sim
#' @title Compute steady-state ARLs of EO-CUSUM control charts using simulation
#' @description Compute steady-state ARLs of EO-CUSUM control charts using simulation.
#'
#' @inheritParams eocusum_arloc_sim
#' @param m Integer. Simulated in-control observations.
#' @param type Character. Default argument is \code{"cond"} for computation of conditional
#' steady-state. Other option is the cyclical steady-state \code{"cycl"}.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template eocusum_ad_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_ad_sim <- function(r, k, h, df, coeff, coeff2, QS = 1, side = "low", type = "cond", m = 50) {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.null(coeff2) || is.na(coeff2) || length(coeff2) != 2) {stop("Model coefficients 'coeff2' must be a numeric vector with two elements")}
  coeff2 <- as.vector(coeff2)
  iside <- switch(as.character(side), low = 1, up = 2)
  if (is.null(iside)) {
    warning("No valid input, using side='low' (deterioration) as default")
    iside <- 1
  }
  QS <- as.numeric(QS)
  if (is.na(QS) || QS <= 0) {stop("QS must a positive numeric value")}
  else if (QS < 1 && iside == 1) {stop("For detecting deterioration (side='low') QS must a positive numeric value >= 1")}
  else if (QS > 1 && iside == 2) {stop("For detecting improvement (side='up') QS must a positive numeric value <= 1")}
  itype <- switch(type, cond = 1, cycl = 2)
  if (is.null(itype)) {
    warning("No valid input, using type='cond' (conditional steady-state) as default")
    itype <- 1
  }
  m <- as.integer(m)
  if (is.na(m) || m < 0) {stop("m must be a positive integer")}
  .eocusum_ad_sim(r, k, h, df, coeff, coeff2, QS, iside, itype, m)
}

#' @name eocusum_arloc_h_sim
#' @title Compute alarm threshold of Out of Control EO-CUSUM control charts using simulation
#' @description Compute alarm threshold (Out of Control ARL) of EO-CUSUM control charts using
#'  simulation.

#' @param L0 Double. Prespecified in-control Average Run Length.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (\code{0/1})
#'  outcome values of each operation.
#' @param k Double. Reference value of the CUSUM control chart. Either \code{0} or a positive
#' value. Can be determined with function \code{\link{optimal_k}}.
#' @param m Integer. Number of simulation runs.
#' @param QS Double. Defines the performance of a surgeon with the odds ratio ratio of death
#' \code{Q}.
#' @param side Character. Default is \code{"low"} to calculate ARL for the upper arm of the V-mask.
#'  If side = \code{"up"}, calculate the lower arm of the V-mask.
#' @param coeff Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model.
#' @param coeff2 Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model of a resampled dataset.
#' @param nc Integer. Number of cores.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given ARL.
#'
#' @details The function \code{eocusum_arloc_h_sim} determines the control limit for given
#' in-control ARL (\code{L0}) by applying a multi-stage search procedure which includes secant rule
#'  and the parallel version of \code{\link{eocusum_arloc_sim}} using \code{\link{mclapply}}.
#'
#' @template eocusum_arloc_h_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_arloc_h_sim <- function(L0, k, df, coeff, coeff2, m = 100, QS = 1, side = "low", nc = 1, jmax = 4, verbose = FALSE) {
  L0 <- as.integer(L0)
  if (is.na(L0) || L0 <= 0) {stop("Given in-control ARL 'L0' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.null(coeff2) || is.na(coeff2) || length(coeff2) != 2) {stop("Model coefficients 'coeff2' must be a numeric vector with two elements")}
  coeff2 <- as.vector(coeff2)
  if (is.null(side)) {
    warning("No valid input, using side='low' (deterioration) as default")
    side <- c("low")
  }
  QS <- as.numeric(QS)
  if (is.na(QS) || QS <= 0) {stop("QS must a positive numeric value")}
  else if (QS < 1 && side == 1) {stop("For detecting deterioration (side='low') QS must a positive numeric value >= 1")}
  else if (QS > 1 && side == 2) {stop("For detecting improvement (side='up') QS must a positive numeric value <= 1")}
  for ( h in 1:10 ) {
    L1 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arloc_sim, h = h, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      L1 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arloc_sim, h = h, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax

  h
}

#' @name eocusum_arl_h_sim
#' @title Compute alarm threshold of EO-CUSUM control charts using simulation
#' @description Compute alarm threshold of EO-CUSUM control charts using simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @param k Double. Reference value of the CUSUM control chart. Either \code{0} or a positive
#' value. Can be determined with function \code{\link{optimal_k}}.
#' @param m Integer. Number of simulation runs.
#' @param side Character. Default is \code{"low"} to calculate ARL for the upper arm of the V-mask.
#'  If side = \code{"up"}, calculate the lower arm of the V-mask.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (0/1)
#' outcome values of each operation.
#' @param coeff Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#' from the binary logistic regression model. For more information see details.
#' @param yemp Logical. If \code{TRUE} use observed outcome value, if \code{FALSE} use estimated
#'  binary logistc regression model.
#' @param nc Integer. Number of cores.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given ARL.
#'
#' @details The function \code{eocusum_arl_h_sim} determines the control limit for given in-control
#'  ARL (\code{L0}) by applying a multi-stage search procedure which includes secant rule and the
#'   parallel version of \code{\link{eocusum_arl_sim}} using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#'
#' @export
#' @examples
#'\donttest{
#'
#'# This function is deprecated. See eocusum_crit_sim() instead.
#'
#'  }
eocusum_arl_h_sim <- function(L0, k, df, coeff, m = 100, yemp = TRUE, side = "low", nc = 1, jmax = 4, verbose = FALSE) {

  .Deprecated("eocusum_arl_h_sim")
  eocusum_crit_sim(L0 = L0, k = k, df = df, coeff = coeff, m = m, yemp = yemp, side = side, nc = nc, jmax = jmax, verbose = verbose)
}

#' @name eocusum_crit_sim
#' @title Compute alarm threshold of EO-CUSUM control charts using simulation
#' @description Compute alarm threshold of EO-CUSUM control charts using simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @param k Double. Reference value of the CUSUM control chart. Either \code{0} or a positive
#' value. Can be determined with function \code{\link{optimal_k}}.
#' @param m Integer. Number of simulation runs.
#' @param side Character. Default is \code{"low"} to calculate ARL for the upper arm of the V-mask.
#'  If side = \code{"up"}, calculate the lower arm of the V-mask.
#' @param df Data Frame. First column are Parsonnet Score values within a range of \code{0} to
#' \code{100} representing the preoperative patient risk. The second column are binary (0/1)
#' outcome values of each operation.
#' @param coeff Numeric Vector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#' from the binary logistic regression model. For more information see details.
#' @param yemp Logical. If \code{TRUE} use observed outcome value, if \code{FALSE} use estimated
#'  binary logistc regression model.
#' @param nc Integer. Number of cores.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given ARL.
#'
#' @template eocusum_crit_sim
#'
#' @details The function \code{eocusum_crit_sim} determines the control limit for given in-control
#'  ARL (\code{L0}) by applying a multi-stage search procedure which includes secant rule and the
#'   parallel version of \code{\link{eocusum_arl_sim}} using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#' @export
eocusum_crit_sim <- function(L0, k, df, coeff, m = 100, yemp = TRUE, side = "low", nc = 1, jmax = 4, verbose = FALSE) {
  L0 <- as.integer(L0)
  if (is.na(L0) || L0 <= 0) {stop("Given in-control ARL 'L0' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.null(side)) {
    warning("No valid input, using side='low' (deterioration) as default")
    side <- c("low")
  }
  for ( h in 1:10 ) {
    L1 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, k = k, h = h, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      L1 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, k = k, h = h, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax

  h
}
