#' @name racusum_discretebeta_arl_sim
#' @title Compute ARLs of RA-CUSUM control charts using simulation
#' @description Compute ARLs of RA-CUSUM control charts using simulation.
#'
#' @param r Integer Vector. Number of runs.
#' @param shape1 Double. Shape parameter \eqn{\alpha}{alpha} \code{> 0} of the beta distribution.
#' @param shape2 Double. Shape parameter \eqn{\beta}{beta} \code{> 0} of the beta distribution.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}.
#' @param coeff Numeric Vector. Estimated intercept and slope coefficients
#'  from a binary logistic regression model.
#' @param RQ Double. Defines the performance of a surgeon with the odds ratio ratio of death.
#' \code{Q}.
#' @param h Double. Control Chart limit for detecting deterioration/improvement.
#' @param rs Integer. Number of intervals between \code{0} and the maximum risk score.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template racusum_discretebeta_arl_sim
#'
#' @author Philipp Wittenberg
#' @export
racusum_discretebeta_arl_sim <- function(r, shape1, shape2, coeff, h, RA = 2, rs = 72, RQ = 1) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_integerish(r, lower = 1, add = arg_checks)
  checkmate::assert_vector(coeff, len = 2, add = arg_checks)
  checkmate::assert_numeric(h, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(rs, lower = 1, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_discretebeta_arl_sim(r, shape1, shape2, coeff, h, RA, rs, RQ)
}

#' @name racusum_discretebeta_crit_sim
#' @title Compute alarm threshold of RA-CUSUM control charts using simulation
#' @description Compute alarm threshold of risk-adjusted cumulative sum control charts using
#'  simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @inheritParams racusum_discretebeta_arl_sim
#' @param nc Integer. Number of cores used for parallel processing. Value is passed to
#' \code{\link{parSapply}}.
#' @param m Integer. Number of simulation runs.
#' @param hmax Integer. Maximum value of \code{h} for the grid search.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given in-control ARL.
#'
#' @template racusum_discretebeta_crit_sim
#'
#' @details Determines the control limit ("\code{h}") for given in-control ARL (\code{"L0"})
#' applying a grid search using \code{\link{racusum_discretebeta_arl_sim}} and
#' \code{\link{parSapply}}.
#'
#' @author Philipp Wittenberg
#'
#' @export
racusum_discretebeta_crit_sim <- function(L0, shape1, shape2, coeff, rs = 72, RA = 2, RQ = 1, nc = 1, hmax = 30, jmax = 4, m = 1e4, verbose = FALSE) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(L0, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(shape1, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(shape2, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_vector(coeff, len = 2, add = arg_checks)
  checkmate::assert_integerish(rs, lower = 1, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(m, lower = 1, add = arg_checks)
  checkmate::assert_integerish(nc, lower = 1, add = arg_checks)
  checkmate::assert_integerish(hmax, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(jmax, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_logical(verbose, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  cl <- parallel::makeCluster(getOption("cl.cores", nc))

  for ( h in 1:hmax ) {
    parallel::clusterExport(cl, c("h", "racusum_discretebeta_arl_sim", "shape1", "shape2", "coeff", "h", "RA", "rs", "RQ", "m"), envir = environment())
    L1 <- mean(parallel::parSapply(cl, 1:m, racusum_discretebeta_arl_sim, shape1 = shape1, shape2 = shape2, coeff = coeff, h = h, RA = RA, rs = rs, RQ = RQ))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      parallel::clusterExport(cl, c("h", "racusum_discretebeta_arl_sim", "shape1", "shape2", "coeff", "h", "RA", "rs", "RQ", "m"), envir = environment())
    L1 <- mean(parallel::parSapply(cl, 1:m, racusum_discretebeta_arl_sim, shape1 = shape1, shape2 = shape2, coeff = coeff, h = h, RA = RA, rs = rs, RQ = RQ))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax
  parallel::stopCluster(cl)
  h
}
