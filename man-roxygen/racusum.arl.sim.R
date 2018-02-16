#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  “Monitoring surgical performance using risk-adjusted cumulative sum charts.”
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441-452.
#'  doi: \href{https://doi.org/10.1093/biostatistics/1.4.441}{10.1093/biostatistics/1.4.441}.
#'
#' @examples
#' \dontrun{
#' require(VLAD2)
#' data(surgeons)
#' S2I <- subset(surgeons, subset=c(surgeon==2 & phase=="I"), select=c(parsonnet, record1))
#' coeff1 <- round(coef(glm(record1~parsonnet, data=S2I, family="binomial")), 3)
#'
#' ## Parallel Simulation 1: y = random (10^4 runs, RA=2)
#' m <- 10^4; h_vec <- 2.7; yemp <- FALSE
#' no_cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(no_cores)
#' parallel::clusterExport(cl, c("h_vec", "racusum.arl.sim", "coeff1", "S2I", "yemp"))
#' time <- system.time( {
#'   ARL <- array(NA, dim=c( length(h_vec), m))
#'   for (h in h_vec) {
#'     ARL[which(h_vec==h), ] <- parallel::parSapply(cl, 1:m, racusum.arl.sim, h=h, coeff=coeff1,
#'                                                  df=S2I, yemp=yemp, USE.NAMES=FALSE) }
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
#' parallel::clusterExport(cl, c("h_vec", "racusum.arl.sim", "coeff1", "S2I"))
#' time <- system.time( {
#'   ARL <- array(NA, dim=c( length(h_vec), m))
#'   for (h in h_vec) {
#'     ARL[which(h_vec==h), ] <- parallel::parSapply(cl, 1:m, racusum.arl.sim, h=h, coeff=coeff1,
#'                                                  df=S2I, USE.NAMES=FALSE) }
#' } )
#' simMean <- apply(ARL, c(1), mean)
#' simSE <- sqrt(apply(ARL, c(1), var)/m)
#' print(list(simMean, simSE, time))
#' parallel::stopCluster(cl)
#' df.sim2 <- data.frame("RA"=2, "h"=h, "ARL"=simMean, "ARLSE"=simSE, "nsim"=m)
#'
#' rbind(df.sim1, df.sim2)
#' }
