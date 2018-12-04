#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  Monitoring surgical performance using risk-adjusted cumulative sum charts.
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441-452.
#'
#' Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455-2473.
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' data("cardiacsurgery", package="spcadjust")
#' df1 <- subset(cardiacsurgery, select=c(Parsonnet, status))
#' coeff1 <- round(coef(glm(status ~ Parsonnet, data=df1, family="binomial")), 3)
#'
#' ## Parallel Simulation 1: y = random (10^4 runs, RA=2)
#' m <- 10^4; h_vec <- 2.7; yemp <- FALSE
#' no_cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(no_cores)
#' parallel::clusterExport(cl, c("h_vec", "racusum_arl_sim", "coeff1", "df1", "yemp"))
#' time <- system.time( {
#'   ARL <- array(NA, dim=c( length(h_vec), m))
#'   for (h in h_vec) {
#'     ARL[which(h_vec==h), ] <- parallel::parSapply(cl, 1:m, racusum_arl_sim, h=h, coeff=coeff1,
#'                                                  df=df1, yemp=yemp, USE.NAMES=FALSE) }
#' } )
#' simMean <- apply(ARL, c(1), mean)
#' simSE <- sqrt(apply(ARL, c(1), var)/m)
#' print(list(simMean, simSE, time))
#' parallel::stopCluster(cl)
#' df.sim1 <- data.frame("RA"=2, "h"=h, "ARL"=simMean, "ARLSE"=simSE, "nsim"=m)
#'
#' ## Parallel Simulation 2: y = empirical (10^4 runs, RA=2)
#' m <- 10^4; h_vec <- 2.7
#' no_cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(no_cores)
#' parallel::clusterExport(cl, c("h_vec", "racusum_arl_sim", "coeff1", "df1"))
#' time <- system.time( {
#'   ARL <- array(NA, dim=c( length(h_vec), m))
#'   for (h in h_vec) {
#'     ARL[which(h_vec==h), ] <- parallel::parSapply(cl, 1:m, racusum_arl_sim, h=h, coeff=coeff1,
#'                                                  df=df1, USE.NAMES=FALSE) }
#' } )
#' simMean <- apply(ARL, c(1), mean)
#' simSE <- sqrt(apply(ARL, c(1), var)/m)
#' print(list(simMean, simSE, time))
#' parallel::stopCluster(cl)
#' df.sim2 <- data.frame("RA"=2, "h"=h, "ARL"=simMean, "ARLSE"=simSE, "nsim"=m)
#'
#' rbind(df.sim1, df.sim2)
#' }
