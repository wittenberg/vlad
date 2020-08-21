#' @name racusum_beta_arl_mc
#' @title Compute ARL of RA-CUSUM control charts assuming patient mix with beta distribution using
#'  Markov chain approximation.
#' @description Compute ARL of risk-adjusted-CUSUM control charts assuming patient mix with beta
#' distribution using Markov chain approximation.
#'
#' @param h Double. \code{h} is the control limit (>\code{0}).
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}. Odds ratio of death under the null hypotheses is \code{1}.
#' @param g0 Double. Estimated intercept coefficient from a binary logistic regression model.
#' @param g1 Double. Estimated slope coefficient from a binary logistic regression model.
#' @param shape1 Double. Shape parameter \eqn{\alpha}{alpha} \code{> 0} of the beta distribution.
#' @param shape2 Double. Shape parameter \eqn{\beta}{beta} \code{> 0} of the beta distribution.
#' @param r Double. Matrix system dimension.
#' @param RQ Double. Defines the performance of a surgeon with the odds ratio ratio of death.
#' @param method Character. If \code{method = "1"} a combination of Sequential Probability Ratio
#'  Test and Toeplitz matrix structure is used to calculate the ARL. \code{"2"} solves a linear
#'  equation system using the classical approach of \emph{Brook and Evans (1972)} to calculate the
#'  ARL.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template racusum_beta_arl_mc
#'
#' @author Philipp Wittenberg
#' @export
racusum_beta_arl_mc <- function(h, RA, g0, g1, shape1, shape2, r = 600, method = 1, RQ = 1) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(h, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(g0, len = 1, add = arg_checks)
  checkmate::assert_numeric(g1, len = 1, add = arg_checks)
  checkmate::assert_numeric(shape1, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(shape2, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(r, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(method, lower = 1, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_beta_arl_mc(h, RA, g0, g1, shape1, shape2, r, method, RQ)
}

#' @name racusum_beta_crit_mc
#' @title Compute alarm threshold of RA-CUSUM control chart assuming patient mix with beta
#' distribution using Markov chain approximation.
#' @description Compute alarm threshold of risk-adjusted CUSUM control chart assuming patient mix
#'  with beta distribution using Markov chain approximation.
#'
#' @param L0 Double. Prespecified Average Run Length.
#' @inheritParams racusum_beta_arl_mc
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{FALSE} a quiet calculation of \code{h} is done. If \code{TRUE}
#'  verbose output of the search procedure (see details) is included.
#'
#' @return Returns a single value which is the control limit \code{h} for a given In-control ARL.
#'
#' @details
#' Determines the control limit for given in-control ARL (\code{"L0"}) using
#' \code{\link{racusum_beta_arl_mc}} by applying a grid search.
#'
#' @template racusum_beta_crit_mc
#'
#' @author Philipp Wittenberg
#' @export
racusum_beta_crit_mc <- function(L0, RA, g0, g1, shape1, shape2, method = 1, r = 600, jmax = 4, verbose = TRUE, RQ = 1) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_integerish(L0, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, lower = 0, add = arg_checks)
  checkmate::assert_numeric(g0, len = 1, add = arg_checks)
  checkmate::assert_numeric(g1, len = 1, add = arg_checks)
  checkmate::assert_numeric(shape1, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(shape2, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(method, lower = 1, add = arg_checks)
  checkmate::assert_numeric(r, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(jmax, lower = 0, add = arg_checks)
  checkmate::assert_logical(verbose, len = 1, add = arg_checks)
  checkmate::assert_numeric(RQ, lower = 0, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_beta_crit_mc(L0, RA, g0, g1, shape1, shape2, method, r, jmax, verbose, RQ)
}
