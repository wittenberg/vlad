#' @name racusum_beta_arl
#' @order 1
#' @title ARL of Beta RA-CUSUM charts
#' @description Compute the ARL of risk-adjusted CUSUM charts assuming a beta distributed
#' patient mix.
#' @param h Double. \code{h} is the control limit (>\code{0}).
#' @param shape1 Double. Shape parameter \eqn{\alpha}{alpha} \code{> 0} of the beta distribution.
#' @param shape2 Double. Shape parameter \eqn{\beta}{beta} \code{> 0} of the beta distribution.
#' @param g0 Double. Estimated intercept coefficient from a binary logistic regression model.
#' @param g1 Double. Estimated slope coefficient from a binary logistic regression model.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}. Odds ratio of death under the null hypotheses is \code{1}.
#' @param RQ Double. Defines the performance of a surgeon with the odds ratio ratio of death.
#' @param r Double. Matrix system dimension.
#' @param method Character. If \code{method = "1"} a combination of Sequential Probability Ratio
#'  Test and Toeplitz matrix structure is used to calculate the ARL. \code{"2"} solves a linear
#'  equation system using the classical approach of \emph{Brook and Evans (1972)} to calculate the
#'  ARL.
#' @return Returns a single value which is the Average Run Length for \code{"racusum_beta_arl_mc"}
#' and \code{"racusum_beta_arl_int"}, and the Run Length for \code{"racusum_beta_arl_sim"}.
#' @references Brook D and Evans DA (1972)
#'  An approach to the probability distribution of CUSUM run length.
#'  \emph{Biometrika}, \strong{59}(3), pp. 539--549
#' @examples
#' \dontrun{
#' library(vlad)
#' ## Markov Chain
#' racusum_beta_arl_mc(h=4.5, shape1=1, shape2=6, g0=-3.6798, g1=0.0768*71, RA=2, r=1e4)
#' ## Full collocation
#' racusum_beta_arl_int(h=4.5, shape1=1, shape2=6, g0=-3.6798, g1=0.0768*71, RA=2, RQ=1, N=150,
#'  pw=FALSE)
#' ## Piece-wise collocation
#' racusum_beta_arl_int(h=4.5, shape1=1, shape2=6, g0=-3.6798, g1=0.0768*71, RA=2, RQ=1, N=49,
#'  pw=TRUE)
#' ## Monte Carlo simulation
#' m <- 1e3
#' RLS <- sapply(1:m, racusum_beta_arl_sim, h=4.5, shape1=1, shape2=6, g0=-3.6798, g1=0.0768,
#' RA = 2, RQ = 1, rs = 71)
#' data.frame(cbind(ARL=mean(RLS), ARLSE=sd(RLS)/sqrt(m)))
#' }
#' @author Philipp Wittenberg
#' @export
racusum_beta_arl_mc <- function(h, shape1, shape2, g0, g1, RA, RQ = 1, r = 600, method = 1) {
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

#' @rdname racusum_beta_arl
#' @order 2
#' @param N Integer. Number of quadrature nodes, dimension of the resulting linear equation system
#'  is equal to \code{N}.
#' @param pw Logical. If \code{FALSE} full collocation is applied. If \code{TRUE} a piece-wise
#'  collocation method is used.
#' @export
racusum_beta_arl_int <- function(h, shape1, shape2, g0, g1, RA, RQ, N, pw) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(h, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(N, lower = 1, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(g0, len = 1, add = arg_checks)
  checkmate::assert_numeric(g1, len = 1, add = arg_checks)
  checkmate::assert_numeric(shape1, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(shape2, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_logical(pw, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_beta_arl_int(h, N, RA, RQ, g0, g1, shape1, shape2, pw)
}

#' @rdname racusum_beta_arl
#' @order 3
#' @param r Integer. Number of runs.
#' @param rs Integer. Maximum risk score.
#' @export
racusum_beta_arl_sim <- function(h, shape1, shape2, g0, g1, r, RA = 2, RQ = 1, rs = 71) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_integerish(r, lower = 1, add = arg_checks)
  checkmate::assert_numeric(g0, len = 1, add = arg_checks)
  checkmate::assert_numeric(g1, len = 1, add = arg_checks)
  coeff <- c(g0, g1)
  checkmate::assert_numeric(h, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(rs, lower = 1, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_beta_arl_sim(r, shape1, shape2, coeff, h, RA, rs, RQ)
}
