#' Variable life-adjusted display (VLAD) and other risk-adjusted control charts
#'
#' The variable life-adjusted display (VLAD) is the first risk-adjusted graphical procedure proposed in the literature
#' for monitoring the performance of a surgeon. It displays the cumulative sum of expected minus observed deaths.
#' Various forms of signaling rules have been developed but they are usually quite complicated.
#' In this package, we establish an equivalence between a VLAD with V-mask and a risk-adjusted cumulative sum (RA-CUSUM) based on the difference between the
#' estimated probability of death and surgical outcome. Average run length analysis based on simulation and markov chain shows that
#' this particular RA-CUSUM chart has similar performance as compared to the established RA-CUSUM chart based on
#' log-likelihood ratio statistic obtained by testing the odds ratio of death. We provide a simple design procedure
#' for determining the V-mask parameters based on a resampling approach. Resampling from a real data set ensures2
#' Finally, we illustrate the monitoring of a real surgeon’s performance using VLAD with V-mask.
#'
#' @docType package
#' @author Philipp Wittenberg <pwitten@hsu-hh.de>
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib vlad
#' @exportPattern "^[[:alpha:]]+"
#' @name vlad-package
NULL
