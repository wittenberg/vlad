#' @name racusum_arl_mc
#' @title Compute ARLs of RA-CUSUM control charts using Markov chain approximation
#' @description Computes the Average Run Length of a risk-adjusted cumulative sum control chart
#' using Markov chain approximation.
#'
#' @param pmix Numeric Matrix. A three column matrix. First column is the risk
#' score distribution. Second column are the predicted probabilities from the risk model. Third
#'  column can be either the predicted probabilities from the risk model or average outcome per
#'  risk score, see examples.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}. Odds ratio of death under the null hypotheses is \code{1}.
#' @param RQ Double. Defines the true performance of a surgeon with the odds ratio ratio of death
#' \code{RQ}. Use \code{RQ = 1} to compute the in-control ARL and other values to compute the
#' out-of-control ARL.
#' @param h Double. \code{h} is the control limit (>\code{0}).
#' @param scaling Double. The \code{scaling} parameter controls the quality of the approximation,
#' larger values achieve higher accuracy but increase the computation burden (larger transition
#' probability matrix).
#' @param rounding Character. If \code{rounding = "p"} a paired rounding implementation of
#' \emph{Knoth et al. (2019)} is used, if \code{rounding = "s"} a simple rounding method of
#' \emph{Steiner et al. (2000)} is used.
#' @param method Character. If \code{method = "Toep"} a combination of Sequential Probability Ratio
#'  Test and Toeplitz matrix structure is used to calculate the ARL. \code{"ToepInv"} computes the
#'  inverted matrix using Toeplitz matrix structure. \code{"BE"} solves a linear equation system
#'  using the classical approach of \emph{Brook and Evans (1972)} to calculate the ARL.
#'
#' @return Returns a single value which is the Average Run Length.
#'
#' @template racusum_arl_mc
#'
#' @author Philipp Wittenberg
#' @export
racusum_arl_mc <- function(pmix, RA, RQ, h, scaling = 600, rounding = "p", method = "Toep") {
  arg_checks <- checkmate::makeAssertCollection()
  pmix <- as.matrix(pmix)
  checkmate::assert_matrix(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(scaling, lower = 1, add = arg_checks)
  rounding <- tolower(rounding)
  checkmate::assert_choice(rounding, c("p", "s"))
  irounding <- switch(rounding, p = 1, s = 2)
  checkmate::assert_choice(method, c("Toep", "ToepInv", "BE"))
  imethod <- switch(method, Toep = 1, ToepInv = 2, BE = 3)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_arl_mc(pmix, RA, RQ, h, scaling, irounding, imethod)
}

#' @name racusum_crit_mc
#' @title Compute alarm threshold of RA-CUSUM control chart using Markov chain approximation
#' @description Computes alarm threshold of a risk-adjusted cumulative sum control chart using
#' Markov chain approximation.
#'
#' @param L0 Double. Prespecified Average Run Length.
#' @param jmax Integer. Number of digits for grid search.
#' @inheritParams racusum_arl_mc
#' @param verbose Logical. If \code{FALSE} a quiet calculation of \code{h} is done. If \code{TRUE}
#'  verbose output of the search procedure is included.
#'
#' @return Returns a single value which is the control limit \code{h} for a given In-control ARL.
#'
#' @details
#' Determines the control limit for given in-control ARL (\code{"L0"}) using
#' \code{\link{racusum_arl_mc}} by applying a grid search.
#'
#' @template racusum_crit_mc
#'
#' @author Philipp Wittenberg
#' @export
racusum_crit_mc <- function(pmix, L0, RA, RQ, scaling = 600, rounding = "p", method = "Toep", jmax = 4, verbose = FALSE) {
  arg_checks <- checkmate::makeAssertCollection()
  pmix <- as.matrix(pmix)
  checkmate::assert_matrix(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_integerish(L0, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(scaling, lower = 1, add = arg_checks)
  rounding <- tolower(rounding)
  checkmate::assert_choice(rounding, c("p", "s"))
  irounding <- switch(rounding, p = 1, s = 2)
  checkmate::assert_choice(method, c("Toep", "ToepInv", "BE"))
  imethod <- switch(method, Toep = 1, ToepInv = 2, BE = 3)
  checkmate::assert_integerish(jmax, lower = 0, add = arg_checks)
  checkmate::assert_logical(verbose, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_crit_mc(pmix, L0, RA, RQ, scaling, irounding, imethod, jmax, verbose)
}
