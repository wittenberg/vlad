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
#' ## Number of simulation runs
#' m <- 10^3
#' ## Number of cores
#' nc <- parallel::detectCores()
#' ## Lower cusum (detecting deterioration)
#' ## k = 0
#' eocusum_arloc_h_sim(L0=370, df=df1, k=0, m=m, side="low", coeff=coeff1, coeff2=coeff2, nc=nc)
#' ## k = kopt
#' QA <- 2
#' # use package function optimal_k to determine k
#' kopt <- optimal_k(QA=QA, parsonnetscores=df1$Parsonnet, coeff=coeff1)
#' eocusum_arloc_h_sim(L0=370, df=df1, k=kopt, m=m, side="low", coeff=coeff1, coeff2=coeff2, nc=nc)
#' ## Upper cusum (detecting improvement)
#' ## k = 0
#' eocusum_arloc_h_sim(L0=370, df=df1, k=0, m=m, side="up", coeff=coeff1, coeff2=coeff2, nc=nc)
#' ## k = kopt
#' QA <- 1/2
#' kopt <- optimal_k(QA=QA, parsonnetscores=df1$Parsonnet, coeff=coeff1)
#' eocusum_arloc_h_sim(L0=370, df=df1, k=kopt, m=m, side="up", coeff=coeff1, coeff2=coeff2, nc=nc)
#' }
