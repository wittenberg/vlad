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
  arg_checks <- checkmate::makeAssertCollection()
  df <- as.data.frame(df)
  checkmate::assert_vector(coeff, len = 2, add = arg_checks)
  checkmate::assert_numeric(R0, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_logical(yemp, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
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
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_integerish(r, len = 1, lower = 1, add = arg_checks)
  # checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(h, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(R0, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .bcusum_arl_sim(r, h, df, R0, RA)
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
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(L0, len = 1, lower = 0, add = arg_checks)
  # checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(R0, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(m, lower = 1, add = arg_checks)
  checkmate::assert_integerish(nc, lower = 1, add = arg_checks)
  checkmate::assert_logical(verbose, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
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

#' @name racusum_ad_sim
#' @title Compute steady-state ARLs of RA-CUSUM control charts using
#' simulation
#' @description Compute steady-state ARLs of risk-adjusted cumulative sum control charts using
#'  simulation.
#' @param r Integer Vector. Number of runs.
#' @param pmix Data Frame. A three column data frame. First column is the operation outcome.
#' Second column are the predicted probabilities from the risk model. Third column can be either the
#'  predicted probabilities from the risk model or average outcome.
#' @param h Double. Control Chart limit for detecting deterioration/improvement.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}. Odds ratio of death under the null hypotheses is \code{1}.
#' @param RQ Double. Defines the true performance of a surgeon with the odds ratio ratio of death
#' \code{RQ}. Use \code{RQ = 1} to compute the in-control ARL and other values to compute the
#' out-of-control ARL.
#' @param m Integer. Simulated in-control observations.
#' @param type Character. Default argument is \code{"cond"} for computation of conditional
#' steady-state. Other option is the cyclical steady-state \code{"cycl"}.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @author Philipp Wittenberg
#'
#' @export
racusum_ad_sim <- function(r, pmix, h, RA = 2, RQ = 1, m = 50, type = "cond") {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_integerish(r, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(h, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(m, lower = 1, add = arg_checks)
  type <- tolower(type)#
  checkmate::assert_choice(type, choices = c("cond", "cycl"), add = arg_checks)
  itype <- switch(type, cond = 1, cycl = 2)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_ad_sim(r, pmix, h, RA, RQ, m, itype)
}

#' @name racusum_arl_sim
#' @title Compute ARLs of RA-CUSUM control charts using simulation
#' @description Computes the Average Run Length of a risk-adjusted cumulative sum control chart
#'  using simulation.
#'
#' @param r Integer Vector. Number of runs.
#' @param pmix Data Frame. A three column data frame. First column is the operation outcome.
#' Second column are the predicted probabilities from the risk model. Third column can be either the
#'  predicted probabilities from the risk model or average outcome.
#' @param h Double. Control Chart limit for detecting deterioration/improvement.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}. Odds ratio of death under the null hypotheses is \code{1}.
#' @param RQ Double. Defines the true performance of a surgeon with the odds ratio ratio of death
#' \code{RQ}. Use \code{RQ = 1} to compute the in-control ARL and other values to compute the
#' out-of-control ARL.
#' @param yemp Logical. If \code{TRUE} use observed outcome value, if \code{FALSE} use estimated
#' binary logistc regression model.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template racusum_arl_sim
#'
#' @author Philipp Wittenberg
#'
#' @export
racusum_arl_sim <- function(r, pmix, h, RA = 2, RQ = 1, yemp = FALSE) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_integerish(r, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(h, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_logical(yemp, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_arl_sim(r, pmix, h, RA, RQ, yemp)
}

#' @name racusum_crit_sim
#' @title Compute alarm threshold of RA-CUSUM control charts using simulation
#' @description Compute alarm threshold of risk-adjusted cumulative sum control charts using
#'  simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @inheritParams racusum_arl_sim
#' @param yemp Logical. If \code{TRUE}, use emirical outcome values, else use model.
#' @param m Integer. Number of simulation runs.
#' @param nc Integer. Number of cores used for parallel processing. Value is passed to
#' \code{\link{parSapply}}.
#' @param hmax Integer. Maximum value of \code{h} for the grid search.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given in-control ARL.
#'
#' @details Determines the control limit ("\code{h}") for given in-control ARL (\code{"L0"})
#' applying a grid search using \code{\link{racusum_arl_sim}} and \code{\link{parSapply}}.
#'
#' @template racusum_crit_sim
#'
#' @author Philipp Wittenberg
#'
#' @export
racusum_crit_sim <- function(L0, pmix, RA = 2, RQ = 1, yemp = FALSE, m = 1e4, nc = 1, hmax = 30, jmax = 4, verbose = FALSE) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(L0, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_logical(yemp, len = 1, add = arg_checks)
  checkmate::assert_integerish(m, lower = 1, add = arg_checks)
  checkmate::assert_integerish(nc, lower = 1, add = arg_checks)
  checkmate::assert_integerish(hmax, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_integerish(jmax, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_logical(verbose, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)

  cl <- parallel::makeCluster(getOption("cl.cores", nc))

  for ( h in 1:hmax ) {
    parallel::clusterExport(cl, c("h", "racusum_arl_sim", "pmix", "RA", "RQ", "yemp", "m"), envir = environment())
    L1 <- mean(parallel::parSapply(cl, 1:m, racusum_arl_sim, pmix = pmix, h = h, RA = RA, RQ = RQ, yemp = yemp))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      parallel::clusterExport(cl, c("h", "racusum_arl_sim", "pmix", "RA", "RQ", "yemp", "m"), envir = environment())
      L1 <- mean(parallel::parSapply(cl, 1:m, racusum_arl_sim, pmix = pmix, h = h, RA = RA, RQ = RQ, yemp = yemp))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax
  parallel::stopCluster(cl)
  h
}

