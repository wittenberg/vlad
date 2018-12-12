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

#' @name optimal_k
#' @title Compute approximate optimal k
#' @description Compute approximate optimal k.
#'
#' @param pmix Data Frame. A three column data frame. First column is the operation outcome.
#' Second column are the predicted probabilities from the risk model. Third
#'  column can be either the predicted probabilities from the risk model or TODO .............
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}. Odds ratio of death under the null hypotheses is \code{1}.
#' @param yemp Logical. If \code{TRUE}, use emirical outcome values, else use model.
#'
#' @return Returns a single value which is the approximate optimal \code{k}.
#'
#' @details Formula deterioration:  \deqn{ k{det} = \frac{R{A} - 1 - log(R{A})}{log(R{A})}\bar{p} , R{A} > 1    }
#'          Formula improvement:    \deqn{ k{imp} = \frac{1 - R{A} + log(R{A})}{log(R{A})}\bar{p} , R{A} < 1    }
#'
#' @template optimal_k
#'
#' @author Philipp Wittenberg
#' @export
optimal_k <- function(pmix, RA, yemp = FALSE) {
  pmix <- as.data.frame(pmix)
  if (is.null(RA) || is.na(RA) || RA <= 0) {stop("QA must be a positive numeric value")}
  RA <- as.numeric(RA)
  if (is.na(yemp) || is.logical(yemp) != "TRUE") {warning("Argument 'yemp' must be logical using TRUE as default value")}
  yemp <- as.logical(yemp)
  .optimal_k(pmix, RA, yemp)
}

#' @name eocusum_arl_sim
#' @title Compute ARLs of EO-CUSUM control charts using simulation
#' @description Compute ARLs of EO-CUSUM control charts using simulation.
#'
#' @param pmix Data Frame. A three column data frame. First column is the operation outcome.
#' Second column are the predicted probabilities from the risk model. Third
#'  column can be either the predicted probabilities from the risk model or TODO .............
#' @param r Integer. Number of of simulation runs.
#' @param k Double. Reference value of the CUSUM control chart. Either \code{0} or a positive
#' value. Can be determined with function \code{\link{optimal_k}}.
#' @param h Double. Decision interval (alarm limit, threshold) of the CUSUM control chart.
#' @param side Character. Default is \code{"low"} to calculate ARL for the upper arm of the V-mask.
#'  If side = \code{"up"}, calculate the lower arm of the V-mask.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template eocusum_arl_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_arl_sim <- function(r, pmix, k, h, RQ = 1, yemp = FALSE, side = "low") {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  pmix <- as.data.frame(pmix)
  if (class(pmix) != "data.frame") {stop("Provide a data frame for argument 'pmix'")}
  else if (ncol(pmix) != 3) {stop("Provide a data frame with three columns for argument 'pmix'")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  RQ <- as.numeric(RQ)
  if (is.na(RQ) || RQ <= 0) {stop("QS must a positive numeric value")}
  else if (RQ < 1 && iside == 1) {stop("For detecting deterioration (side='low') RQ must a positive numeric value >= 1")}
  else if (RQ > 1 && iside == 2) {stop("For detecting improvement (side='up') RQ must a positive numeric value <= 1")}
  if (is.na(yemp) || is.logical(yemp) != "TRUE") {warning("Argument 'yemp' must be logical using TRUE as default value")}
  yemp <- as.logical(yemp)
  iside <- switch(as.character(side), low = 1, up = 2)
  if (is.null(iside)) {
    warning("No valid input, using side='low' (deterioration) as default")
    iside <- 1
  }
  .eocusum_arl_sim(r, pmix, k, h, RQ, yemp, iside)
}

#' @name eocusum_ad_sim
#' @title Compute steady-state ARLs of EO-CUSUM control charts using simulation
#' @description Compute steady-state ARLs of EO-CUSUM control charts using simulation.
#'
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
eocusum_ad_sim <- function(r, k, h, pmix, RQ = 1, side = "low", type = "cond", m = 50) {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  pmix <- as.data.frame(pmix)
  if (class(pmix) != "data.frame") {stop("Provide a data frame for argument 'pmix'")}
  else if (ncol(pmix) != 3) {stop("Provide a data frame with three columns for argument 'pmix'")}
  iside <- switch(as.character(side), low = 1, up = 2)
  if (is.null(iside)) {
    warning("No valid input, using side='low' (deterioration) as default")
    iside <- 1
  }
  RQ <- as.numeric(RQ)
  if (is.na(RQ) || RQ <= 0) {stop("QS must a positive numeric value")}
  else if (RQ < 1 && iside == 1) {stop("For detecting deterioration (side='low') QS must a positive numeric value >= 1")}
  else if (RQ > 1 && iside == 2) {stop("For detecting improvement (side='up') QS must a positive numeric value <= 1")}
  itype <- switch(type, cond = 1, cycl = 2)
  if (is.null(itype)) {
    warning("No valid input, using type='cond' (conditional steady-state) as default")
    itype <- 1
  }
  m <- as.integer(m)
  if (is.na(m) || m < 0) {stop("m must be a positive integer")}
  .eocusum_ad_sim(r, pmix, k, h, RQ, iside, itype, m)
}

#' @name eocusum_crit_sim
#' @title Compute alarm threshold of EO-CUSUM control charts using simulation
#' @description Compute alarm threshold of EO-CUSUM control charts using simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @inheritParams eocusum_arl_sim
#' @param yemp Logical. If \code{TRUE}, use emirical outcome values, else use model.
#' @param m Integer. Number of simulation runs.
#' @param nc Integer. Number of cores used for parallel processing.
#' @param hmax Integer. Maximum value of \code{h} for the grid search.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given in-control ARL.
#'
#' @details Determines the control limit ("\code{h}") for given in-control ARL (\code{"L0"})
#' applying a grid search using \code{\link{eocusum_arl_sim}} and \code{\link{parSapply}}.
#'
#'
#' @author Philipp Wittenberg
#'
#' @export
eocusum_crit_sim <- function(L0, pmix, k, RQ = 1, side = "low", yemp = FALSE, m = 1e4, nc = 1, hmax = 30, jmax = 4, verbose = FALSE) {
  L0 <- as.integer(L0)
  if (is.na(L0) || L0 <= 0) {stop("Given in-control ARL 'L0' must be a positive integer")}
  pmix <- as.data.frame(pmix)
  if (class(pmix) != "data.frame") {stop("Provide a data frame for argument 'pmix'")}
  else if (ncol(pmix) != 3) {stop("Provide a data frame with three columns for argument 'pmix'")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  RQ <- as.numeric(RQ)
  if (is.na(RQ) || RQ <= 0) {stop("QS must be a positive numeric value")}
  else if (RQ < 1 && iside == 1) {stop("For detecting deterioration (side='low') QS must a positive numeric value >= 1")}
  else if (RQ > 1 && iside == 2) {stop("For detecting improvement (side='up') QS must a positive numeric value <= 1")}
  iside <- switch(as.character(side), low = 1, up = 2)
  if (is.null(iside)) {
    warning("No valid input, using side='low' (deterioration) as default")
    iside <- 1
  }
  if (is.na(yemp) || is.logical(yemp) != "TRUE") {warning("Argument 'yemp' must be logical using TRUE as default value")}
  yemp <- as.logical(yemp)
  m <- as.integer(m)
  if (is.na(m) || m < 0) {stop("m must be a positive integer")}
  nc <- as.integer(nc)
  if (is.na(nc) || nc <= 0) {stop("Number of cores 'nc' must be a positive integer")}
  hmax <- as.integer(hmax)
  if (is.na(hmax) || hmax <= 0) {stop("Maximum value for grid search 'hmax' must be a positive integer")}
  jmax <- as.integer(jmax)
  if (is.na(jmax) || jmax <= 0) {stop("Number of digits for grid search 'jmax' must be a positive integer")}

  verbose <- as.logical(verbose)

  cl <- parallel::makeCluster(getOption("cl.cores", nc))

  for ( h in 1:hmax ) {
    parallel::clusterExport(cl, c("h", "eocusum_arl_sim", "pmix", "k", "RQ", "iside", "yemp"), envir = environment())
    L1 <- mean(parallel::parSapply(cl, 1:m, eocusum_arl_sim, pmix = pmix, k = k, h = h, RQ = RQ, side = iside, yemp = yemp))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      parallel::clusterExport(cl, c("h", "eocusum_arl_sim", "pmix", "k", "RQ", "iside", "yemp"), envir = environment())
      L1 <- mean(parallel::parSapply(cl, 1:m, eocusum_arl_sim, pmix = pmix, k = k, h = h, RQ = RQ, side = iside, yemp = yemp))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax
  parallel::stopCluster(cl)
  h
}
