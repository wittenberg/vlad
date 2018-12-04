#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  Monitoring surgical performance using risk-adjusted cumulative sum charts.
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441--452.
#'
#' Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455--2473.
#'
#' @examples
#' \dontrun{
#' library("vlad"); library("ggplot2")
#' ## Set seed for reproducibility
#' RNGkind("L'Ecuyer-CMRG")
#' ## Datasets
#' data("cardiacsurgery", package = "spcadjust")
#' s5000 <- dplyr::sample_n(cardiacsurgery, size = 5000, replace = TRUE)
#' df1 <- subset(cardiacsurgery, select = c(Parsonnet, status))
#' df2 <- subset(s5000, select = c(Parsonnet, status))
#'
#' ## Estimate coefficients from logit model
#' coeff1 <- round(coef(glm(status ~ Parsonnet, data = df1, family = "binomial")), 3)
#' coeff2 <- round(coef(glm(status ~ Parsonnet, data = df2, family = "binomial")), 3)
#'
#' ## Number of simulation runs
#' m <- 10^3
#'
#' ## Deterioration RA=2:
#' ## 1. Determine critical value for given ARL
#' h0 <- racusum_arloc_h_sim(L0 = 370, df = df1, coeff = coeff1, coeff2 = coeff2, m = m, RA = 2,
#' nc = 6)
#' ## 2. Compute Out of Control ARL
#' RQ <- seq(1, 4, 0.1)
#' rl <- array(NA, dim = c(m, length(RQ)))
#' RLS <- sapply(RQ, function(i) {
#'   cat("RQ: ", i, "\n" )
#'   rl[, i] <- do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h0, df = df1, RA = 2,
#'   RQ = i, coeff = coeff1, coeff2 = coeff2, mc.cores = 6))
#' })
#' df3 <- data.frame(cbind(RQ, "ARL" = apply(RLS, 2, mean), "ARLSE" = apply(RLS, 2, mean)/sqrt(m)))
#' ggplot(df3, aes(RQ, ARL)) + geom_line() + theme_classic()
#'
#' ## Improvement RA=1/2:
#' ## 1. Determine critical value for given ARL
#' h0 <- racusum_arloc_h_sim(L0 = 370, df = df1, coeff = coeff1, coeff2 = coeff2, m = m, RA = 1/2,
#'                           nc = 6)
#' ## 2. Compute Out of Control ARL
#' RQ <- seq(1/4, 1, 1/40)
#' rl <- array(NA, dim = c(m, length(RQ)))
#' RLS <- sapply(RQ, function(i) {
#'   cat("RQ: ", i, "\n" )
#'   rl[, i] <- do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h = h0, df = df1, RA = 1/2,
#'                                            RQ = i, coeff = coeff1, coeff2 = coeff2,
#'                                            mc.cores = 6))
#' })
#' df4 <- data.frame(cbind(RQ, "ARL" = apply(RLS, 2, mean), "ARLSE" = apply(RLS, 2, mean)/sqrt(m)))
#' ggplot(df4, aes(RQ, ARL)) + geom_line() + theme_classic()
#' }
