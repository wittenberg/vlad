#' @name racusum_beta_crit
#' @title Alarm thresholds of Beta RA-CUSUM charts
#' @description Compute alarm threshold of risk-adjusted CUSUM charts assuming a beta distributed
#' patient mix.
#' @param L0 Double. Prespecified Average Run Length.
#' @param shape1 Double. Shape parameter \eqn{\alpha}{alpha} \code{> 0} of the beta distribution.
#' @param shape2 Double. Shape parameter \eqn{\beta}{beta} \code{> 0} of the beta distribution.
#' @param g0 Double. Estimated intercept coefficient from a binary logistic regression model.
#' @param g1 Double. Estimated slope coefficient from a binary logistic regression model.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}. Odds ratio of death under the null hypotheses is \code{1}.
#' @param r Double. Matrix system dimension.
#' @param RQ Double. Defines the true performance of a surgeon with the odds ratio ratio of death
#' \code{RQ}. Use \code{RQ = 1} to compute the in-control ARL and other values to compute the
#' out-of-control ARL.
#' @param method Character. If \code{method = "1"} a combination of Sequential Probability Ratio
#'  Test and Toeplitz matrix structure is used to calculate the ARL. \code{"2"} solves a linear
#'  equation system using the classical approach of \emph{Brook and Evans (1972)} to calculate the
#'  ARL.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{FALSE} a quiet calculation of \code{h} is done. If \code{TRUE}
#'  verbose output of the search procedure (see details) is included.
#' @return Returns a single value which is the control limit \code{h} for a given In-control ARL.
#' @details Determines the control limit ("\code{h}") for a given in-control ARL (\code{"L0"}) using
#' \code{\link{racusum_beta_arl_mc}} or \code{\link{racusum_beta_arl_sim}} and
#' \code{\link{parSapply}} by applying a grid search.
#' @references Brook D and Evans DA (1972)
#'  An approach to the probability distribution of CUSUM run length.
#'  \emph{Biometrika}, \strong{59}(3), pp. 539--549
#' @examples
#' \dontrun{
#' library(vlad)
#' ## Markov Chain
#' racusum_beta_crit_mc(L0=7500, shape1=.61, shape2=4.09, g0=-3.6798, g1=0.0768*71, RA=2, RQ=1,
#'  r=1e3)
#' ## Monte Carlo simulation
#' racusum_beta_crit_sim(L0=7500, shape1=.61, shape2=4.09, g0=-3.6798, g1=0.0768, RA = 2, RQ = 1,
#' rs = 71, verbose=TRUE, m=1e3)
#' }
#' @export
racusum_beta_crit_mc <- function(L0, shape1, shape2, g0, g1, RA, RQ = 1, method = 1, r = 600, jmax = 4, verbose = TRUE) {
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

#' @rdname racusum_beta_crit
#' @param nc Integer. Number of cores used for parallel processing. Value is passed to
#' \code{\link{parSapply}}.
#' @param m Integer. Number of simulation runs.
#' @param rs Integer. Maximum risk score.
#' @param hmax Integer. Maximum value of \code{h} for the grid search.
#' @export
racusum_beta_crit_sim <- function(L0, shape1, shape2, g0, g1, RA = 2, RQ = 1, nc = 1, rs = 71, hmax = 30, jmax = 4, m = 1e4, verbose = FALSE) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(L0, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(shape1, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(shape2, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(g0, len = 1, add = arg_checks)
  checkmate::assert_numeric(g1, len = 1, add = arg_checks)
  coeff <- c(g0, g1)
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
    parallel::clusterExport(cl, c("h", "racusum_beta_arl_sim", "shape1", "shape2", "g0", "g1", "h", "RA", "rs", "RQ", "m"), envir = environment())
    L1 <- mean(parallel::parSapply(cl, 1:m, racusum_beta_arl_sim, shape1 = shape1, shape2 = shape2, g0 = g0, g1 = g1, h = h, RA = RA, rs = rs, RQ = RQ))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      parallel::clusterExport(cl, c("h", "racusum_beta_arl_sim", "shape1", "shape2", "g0", "g1", "h", "RA", "rs", "RQ", "m"), envir = environment())
    L1 <- mean(parallel::parSapply(cl, 1:m, racusum_beta_arl_sim, shape1 = shape1, shape2 = shape2, g0 = g0, g1 = g1, h = h, RA = RA, rs = rs, RQ = RQ))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax
  parallel::stopCluster(cl)
  h
}
