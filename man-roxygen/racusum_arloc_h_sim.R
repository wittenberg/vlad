#' @examples
#' \dontrun{
#' require("vlad")
#' # Set seed for reproducibility
#' RNGkind("L'Ecuyer-CMRG")
#' set.seed(1234)
#' parallel::mc.reset.stream()
#' # Datasets
#' data("cardiacsurgery")
#' s5000 <- dplyr::sample_n(cardiacsurgery, size=5000, replace=TRUE)
#' df1 <- subset(cardiacsurgery, select=c(Parsonnet, status))
#' df2 <- subset(s5000, select=c(Parsonnet, status))
#'
#' # Estimate coefficients from logit model
#' coeff1 <- round(coef(glm(status~Parsonnet, data=df1, family="binomial")), 3)
#' coeff2 <- round(coef(glm(status~Parsonnet, data=df2, family="binomial")), 3)
#'
#' # Number of simulation runs
#' m <- 10^3
#'
#' # Deterioration:
#' # 1. Determine critical value for given ARL
#' racusum_arloc_h_sim(L0=370, df=df1, coeff=coeff1, coeff2=coeff2, m=m, RA=2, nc=6)
#' # h=2.030933
#'
#' # 2. Determine ARL and Standard Error
#' RLS <- do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h=2.035, df=df1, RA=2, coeff=coeff1,
#'                                      coeff2=coeff2, mc.cores=6))
#' data.frame(cbind("ARL"=mean(RLS), "ARLSE"=sd(RLS)/sqrt(m)))
#' # ARL=371.125; ARLSE=11.36053
#'
#' # Improvement:
#' # 1. Determine critical value for given ARL
#' racusum_arloc_h_sim(L0=370, df=df1, coeff=coeff1, coeff2=coeff2, m=m, RA=1/2, nc=6)
#' # h=1.710999
#' #
#' # 2. Determine ARL and Standard Error
#' RLS <- do.call(c, parallel::mclapply(1:m, racusum_arloc_sim, h=1.760, df=df1, RA=1/2, coeff=coeff1,
#'                                      coeff2=coeff2, mc.cores=6))
#' data.frame(cbind("ARL"=mean(RLS), "ARLSE"=sd(RLS)/sqrt(m)))
#' # ARL=399.613; ARLSE=10.7601
#' }
