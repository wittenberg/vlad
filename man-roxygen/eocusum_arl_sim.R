#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  “Monitoring surgical performance using risk-adjusted cumulative sum charts.”
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441-452.
#'  doi: \href{https://doi.org/10.1093/biostatistics/1.4.441}{10.1093/biostatistics/1.4.441}.
#'
#' @examples
#' \dontrun{
#' require("vlad")
#' require("spcadjust")
#' # Set seed for reproducibility
#' RNGkind("L'Ecuyer-CMRG")
#' # Datasets
#' data("cardiacsurgery")
#' df1 <- subset(cardiacsurgery, select=c(Parsonnet, status))
#'
#' # Estimate coefficients from logit model
#' coeff1 <- round(coef(glm(status~Parsonnet, data=df1, family="binomial")), 3)
#'
#' # Number of cores
#'
#' nc <- parallel::detectCores()
#'
#' # Number of simulation runs
#' m <- 10^3
#' # Deterioration:
#' kopt <- optimal_k(QA=2, parsonnetscores=df1$Parsonnet, coeff=coeff1)
#' # 1. Determine critical value for given ARL
#' hdet <- eocusum_arl_h_sim(L0=370, df=df1, k=kopt, coeff=coeff1, m=m, nc=nc, side="low")
#' 2. Determine ARL and Standard Error
#' RLS <- do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, h=hdet, k=kopt, df=df1, side="low",
#'                                      coeff=coeff1, mc.cores=nc))
#' data.frame(cbind("ARL"=mean(RLS), "ARLSE"=sd(RLS)/sqrt(m)))
#'
#' # Improvement:
#' kopt <- optimal_k(QA=1/2, parsonnetscores=df1$Parsonnet, coeff=coeff1)
#' # 1. Determine critical value for given ARL
#' hup <- eocusum_arl_h_sim(L0=370, df=df1, k=kopt, coeff=coeff1, m=m, nc=nc, side="up")
#' #
#' # 2. Determine ARL and Standard Error
#' RLS <- do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, h=hup, k=kopt, df=df1, side="up",
#'                                      coeff=coeff1, mc.cores=nc))
#' data.frame(cbind("ARL"=mean(RLS), "ARLSE"=sd(RLS)/sqrt(m)))
#' }
