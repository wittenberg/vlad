#' @references Wittenberg P, Gan FF, Knoth S (2018).
#' “A simple signaling rule for variable life‐adjusted display derived from
#' an equivalent risk‐adjusted CUSUM chart.”
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455-2473.
#' doi: \doi{10.1002/sim.7647}.
#'
#' @examples
#' \dontrun{
#' ## Datasets
#' data("cardiacsurgery", package = "spcadjust")
#' s5000 <- dplyr::sample_n(cardiacsurgery, size = 5000, replace = TRUE)
#' df1 <- subset(cardiacsurgery, select=c(Parsonnet, status))
#' df2 <- subset(s5000, select=c(Parsonnet, status))
#' ## estimate coefficients from logit model
#' coeff1 <- round(coef(glm(status~Parsonnet, data=df1, family="binomial")), 3)
#' coeff2 <- round(coef(glm(status~Parsonnet, data=df2, family="binomial")), 3)
#'
#' ## Serial simulation
#' ## set seed for reproducibility
#' RNGkind("L'Ecuyer-CMRG")
#' m <- 10^3
#' RLS <- do.call(c, lapply(1:m, eocusum_arloc_sim, h=4.498, k=0, df=df1, side="low", coeff=coeff1,
#'                          coeff2=coeff2))
#' data.frame(cbind("ARL"=mean(RLS), "ARLSE"=sd(RLS)/sqrt(m)))
#' ## ARL=366.697; ARLSE=9.457748
#' ## Parallel simulation (FORK)
#' ## set seed for reproducibility
#' RNGkind("L'Ecuyer-CMRG")
#' RLS <- simplify2array(parallel::mclapply(1:m, eocusum_arloc_sim, h=4.498, k=0, df=df1, side="low",
#'                                          coeff=coeff1, coeff2=coeff2,
#'                                          mc.cores=parallel::detectCores()))
#' data.frame(cbind("ARL"=mean(RLS), "ARLSE"=sd(RLS)/sqrt(m)))
#'
#' ## Parallel simulation (PSOCK)
#' ## set seed for reproducibility
#' RNGkind("L'Ecuyer-CMRG")
#' no_cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(no_cores)
#' side <- 1
#' h_vec <- 4.498
#' QS_vec <- 1
#' m <- 10^3
#' k <- 0
#' parallel::clusterExport(cl, c("h_vec", "eocusum_arloc_sim", "df1", "coeff1", "coeff2",
#'                               "QS_vec", "side", "k"))
#' time <- system.time( {
#' RLS <- array(NA, dim=c( length(QS_vec), length(h_vec), m))
#' for (h in h_vec) {
#'   for (QS in QS_vec) {
#'     cat(h, " ", QS, "\n")
#'     RLS[which(QS_vec==QS), which(h==h_vec), ] <- parallel::parSapply(cl, 1:m, eocusum_arloc_sim,
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
#' }
