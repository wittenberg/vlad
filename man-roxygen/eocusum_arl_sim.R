#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  “Monitoring surgical performance using risk-adjusted cumulative sum charts.”
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441-452.
#'  doi: \href{https://doi.org/10.1093/biostatistics/1.4.441}{10.1093/biostatistics/1.4.441}.
#'
#' @examples
#' \dontrun{
#' library("vlad")
#' library("spcadjust")
#' ## Datasets
#' data("cardiacsurgery")
#' s5000 <- dplyr::sample_n(cardiacsurgery, size=5000, replace=TRUE)
#' df1 <- subset(cardiacsurgery, select=c(Parsonnet, status))
#' df2 <- subset(s5000, select=c(Parsonnet, status))
#' ## estimate coefficients from logit model
#' coeff1 <- round(coef(glm(status~Parsonnet, data=df1, family="binomial")), 3)
#' coeff2 <- round(coef(glm(status~Parsonnet, data=df2, family="binomial")), 3)
#'
#' ## Serial simulation
#' RNGkind("L'Ecuyer-CMRG")
#' m <- 10^3
#' kopt <- optimal_k(QA=2, parsonnetscores=df1$Parsonnet, coeff=coeff1)
#' #eocusum_arloc_h_sim(L0=370, df=df1, k=kopt, m=m, side="low", coeff=coeff1, coeff2=coeff2, nc=nc)
#' RLS <- do.call(c, lapply(1:m, eocusum_arloc_sim, h=2.626, k=kopt, df=df1, side="low", coeff=coeff1,
#'                          coeff2=coeff2))
#' data.frame(cbind(ARL=mean(RLS), ARLSE=sd(RLS)/sqrt(m)))
#' ## Parallel simulation (FORK)
#' RNGkind("L'Ecuyer-CMRG")
#' m <- 10^3
#' kopt <- optimal_k(QA=2, parsonnetscores=df1$Parsonnet, coeff=coeff1)
#' RLS <- simplify2array(parallel::mclapply(1:m, eocusum_arloc_sim, h=2.626, k=kopt, df=df1, side="low",
#'                                          coeff=coeff1, coeff2=coeff2,
#'                                          mc.cores=parallel::detectCores()))
#' data.frame(cbind(ARL=mean(RLS), ARLSE=sd(RLS)/sqrt(m)))
#'
#' ## Parallel simulation (PSOCK)
#' RNGkind("L'Ecuyer-CMRG")
#' no_cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(no_cores)
#' side <- "low"
#' h_vec <- 2.626
#' QS_vec <- 1
#' m <- 10^3
#' k <- optimal_k(QA=2, parsonnetscores=df1$Parsonnet, coeff=coeff1)
#' parallel::clusterExport(cl, c("h_vec", "eocusum_arloc_sim", "df1", "coeff1", "coeff2",
#'                               "QS_vec", "side", "k"))
#' time <- system.time( {
#'   RLS <- array(NA, dim=c( length(QS_vec), length(h_vec), m))
#'   for (h in h_vec) {
#'     for (QS in QS_vec) {
#'       cat(h, " ", QS, "\n")
#'       RLS[which(QS_vec==QS), which(h==h_vec), ] <- parallel::parSapply(cl, 1:m, eocusum_arloc_sim,
#'                                                                        side=side, QS=QS, h=h, k=k,
#'                                                                        df=df1,  coeff=coeff1,
#'                                                                        coeff2=coeff2,
#'                                                                        USE.NAMES=FALSE)
#'     }
#'   }
#' } )
#' ARL <- apply(RLS, c(1, 2), mean)
#' ARLSE <- sqrt(apply(RLS, c(1, 2), var)/m)
#' print(list(ARL, ARLSE, time))
#' parallel::stopCluster(cl)
#' }
