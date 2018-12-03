#' Variable Life Adjusted Display and Other Risk-Adjusted Quality Control Charts
#'
#' Contains functions to set up risk-adjusted quality control charts
#' in health care. For the variable life adjusted display (VLAD) proposed by
#' Lovegrove et al. (1997) <doi:10.1016/S0140-6736(97)06507-0> signaling rules
#' derived in Wittenberg et al. (2018) <doi: 10.1002/sim.7647> are implemented.
#' Additionally, for the risk-adjusted cumulative sum chart based on log-likelihood
#' ratio statistic introduced by Steiner et al. (2000) <doi:10.1093/biostatistics/1.4.441>
#'   average run length and control limits can be computed.
#' @docType package
#' @importFrom Rcpp evalCpp
#' @useDynLib vlad, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#' @name vlad-package
NULL
