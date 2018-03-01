#' Variable Life Adjusted Display
#'
#' Contains functions to set up risk-adjusted quality control charts
#' in health care. For the variable life adjusted display (VLAD) proposed by
#' Lovegrove et al. (1997) <doi:10.1016/S0140-6736(97)06507-0> and the
#' risk-adjusted cumulative sum chart based on log-likelihood ratio statistic
#' introduced by Steiner et al. (2000) <doi:10.1093/biostatistics/1.4.441> the
#' average run length and control limits can be computed.
#' @docType package
#' @author Philipp Wittenberg and Sven Knoth
#' @importFrom Rcpp evalCpp
#' @useDynLib vlad, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#' @name vlad-package
NULL
