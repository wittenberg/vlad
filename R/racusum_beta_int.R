#' @name racusum_beta_arl_int
#' @title Compute ARL of RA-CUSUM control charts assuming a beta distributed patient mix with using
#'  collocation methods.
#' @description Compute ARL of RA-CUSUM control charts assuming a beta distributed patient mix with using
#'  collocation methods.
#'
#' @param h Double. \code{h} is the control limit (>\code{0}).
#' @param N Integer. Number of quadrature nodes, dimension of the resulting linear equation system
#'  is equal TODO ???.
#' @param RA ...
#' @param RQ ...
#' @param g0 Double. ...
#' @param g1 Double. ...
#' @param shape1 Double. Shape parameter \eqn{\alpha}{alpha} \code{> 0} of the beta distribution.
#' @param shape2 Double. Shape parameter \eqn{\beta}{beta} \code{> 0} of the beta distribution.
#' @param pw Logical. Pice-wise collocation yes or no.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template racusum_beta_arl_int
#'
#' @author Philipp Wittenberg
#' @export
racusum_beta_arl_int <- function(h, N, RA, RQ, g0, g1, shape1, shape2, pw) {
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
