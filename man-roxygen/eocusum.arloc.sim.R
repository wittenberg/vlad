#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  “Monitoring surgical performance using risk-adjusted cumulative sum charts.”
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441-452.
#'  doi: \href{https://doi.org/10.1093/biostatistics/1.4.441}{10.1093/biostatistics/1.4.441}.
#'
#' @examples
#' \dontrun{
#' rm(list=ls())
#' require(VLAD2)
#' data("surgeons"); data("s5000")
#' df1 <- subset(surgeons, select=c(parsonnet, record1))
#' df2 <- subset(s5000, select=c(parsonnet, record1))
#' ## estimate coefficients from logit model
#' coeff1 <- round(coef(glm(record1~parsonnet, data=df1, family="binomial")), 3)
#' coeff2 <- round(coef(glm(record1~parsonnet, data=df2, family="binomial")), 3)
#'
#' ## Serial simulation
#' ## set seed for reproducibility
#' RNGkind("L'Ecuyer-CMRG") # RNGkind("Mersenne-Twister")
#' set.seed(1234)
#' m <- 10^3
#' RLS <- do.call(c, lapply(1:m, eocusum.arloc.sim, h=4.498, k=0, df=df1, side=1, coeff=coeff1,
#'                          coeff2=coeff2))
#' data.frame(cbind("ARL"=mean(RLS), "ARLSE"=sd(RLS)/sqrt(m)))
#' ## ARL=366.697; ARLSE=9.457748
#' ## Parallel simulation (FORK)
#' ## set seed for reproducibility
#' RNGkind("L'Ecuyer-CMRG")
#' set.seed(1234); parallel::mc.reset.stream()
#' m <- 10^3
#' RLS <- simplify2array(parallel::mclapply(1:m, eocusum.arloc.sim, h=4.498, k=0, df=df1, side=1,
#'                                          coeff=coeff1, coeff2=coeff2,
#'                                          mc.cores=parallel::detectCores()))
#' data.frame(cbind("ARL"=mean(RLS), "ARLSE"=sd(RLS)/sqrt(m)))
#' ## ARL=372.133; ARLSE=9.967791
#'
#' ## Parallel simulation (PSOCK)
#' ## set seed for reproducibility
#' RNGkind("L'Ecuyer-CMRG")
#' no_cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(no_cores)
#' parallel::clusterSetRNGStream(cl, iseed = 1234)
#' side <- 1
#' h_vec <- 4.498
#' QS_vec <- 1
#' m <- 10^3
#' k <- 0
#' parallel::clusterExport(cl, c("h_vec", "eocusum.arloc.sim", "df1", "coeff1", "coeff2",
#'                               "QS_vec", "side", "k"))
#' time <- system.time( {
#' RLS <- array(NA, dim=c( length(QS_vec), length(h_vec), m))
#' for (h in h_vec) {
#'   for (QS in QS_vec) {
#'     cat(h, " ", QS, "\n")
#'     RLS[which(QS_vec==QS), which(h==h_vec), ] <- parallel::parSapply(cl, 1:m, eocusum.arloc.sim,
#'                                                                      side=side, QS=QS, h=h, k=k,
#'                                                                      df=df1,  coeff=coeff1,
#'                                                                      coeff2=coeff2,
#'                                                                       USE.NAMES=FALSE)
#'     }
#'   }
#' } )
#' ARL <- apply(RLS, c(1, 2), mean)
#' ARLSE <- sqrt(apply(RLS, c(1, 2), var)/m)
#' print(list(ARL, ARLSE, time))
#' parallel::stopCluster(cl)
#' ## ARL=371.361; ARLSE=9.899058
#' }
