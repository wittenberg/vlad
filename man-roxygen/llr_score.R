#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  Monitoring surgical performance using risk-adjusted cumulative sum charts.
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441--452.
#'
#' Steiner S (2014). Risk-Adjusted Monitoring of Outcomes in Health Care.
#' In Lawless JF (ed.), \emph{Statistics in Action}, pp. 225--242. Informa UK Limited.
#'
#' Rigdon SE and Fricker RD (2015). Health Surveillance.
#' In Chen DG and Wilson J (eds) \emph{Innovative Statistical Methods for Public Health Data},
#' pp. 203--249. Springer, Cham.
#'
#' @examples
#' ## see Steiner et al. (2000) p. 446 or Steiner (2014) p. 234
#' coeff <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
#' ## Log-likelihood ratio scores for detecting an increase in the failure rate:
#' ## low risk patients with a Parsonnet score of zero
#'
#' llr_score(df = data.frame(as.integer(0), 0), coeff = coeff, RA = 2)
#' llr_score(df = data.frame(as.integer(0), 1), coeff = coeff, RA = 2)
#'
#' ## higher risk patients with a Parsonnet score of 50
#' llr_score(df = data.frame(as.integer(50), 0), coeff = coeff, RA = 2)
#' llr_score(df = data.frame(as.integer(50), 1), coeff = coeff, RA = 2)
#'
#' ## see Steiner (2014) p. 234
#' ## Log-likelihood ratio scores for detecting an decrease in the failure rate:
#' ## low risk patients with a Parsonnet score of zero
#' llr_score(df = data.frame(as.integer(0), 0), coeff = coeff, RA = 1/2)
#' llr_score(df = data.frame(as.integer(0), 1), coeff = coeff, RA = 1/2)
#'
#' ## higher risk patients with a Parsonnet score of 50
#' llr_score(df = data.frame(as.integer(50), 0), coeff = coeff, RA = 1/2)
#' llr_score(df = data.frame(as.integer(50), 1), coeff = coeff, RA = 1/2)
#'
#' ## see Rigdon and Fricker p. 225 and 226
#' ## detecting an increase in the failure rate:
#' coeff <- c("(Intercept)" = -3.67, "Parsonnet" = 0.077)
#' df <- data.frame(Parsonnet = c(19L, 19L, 0L, 0L), status = c(0, 1, 0, 1))
#' lapply(seq_along(df$Parsonnet), function(i) round(llr_score(df = df[i, ], coeff = coeff,
#'  RA = 2), 4))
#'
#' ## detecting an decrease in the failure rate:
#' round(llr_score(df = data.frame(19L, 0), coeff = coeff, RA = 1/2), 5)
#'
