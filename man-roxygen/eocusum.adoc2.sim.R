#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#' “Monitoring surgical performance using risk-adjusted cumulative sum charts.”
#' \emph{Biostatistics}, \strong{1}(4), pp. 441-452.
#' doi: \href{https://doi.org/10.1093/biostatistics/1.4.441}{10.1093/biostatistics/1.4.441}.
#'
#' Taylor HM (1968). “The Economic Design of Cumulative Sum Control Charts.”
#' \emph{Technometrics}, \strong{10}(3), pp. 479-488.
#' doi: \href{https://doi.org/10.1080/00401706.1968.10490595}{10.1080/00401706.1968.10490595}.
#'
#' Crosier R (1986). “A new two-sided cumulative quality control scheme.”
#' \emph{Technometrics}, \strong{28}(3), pp. 187-194.
#' doi: \href{https://doi.org/10.2307/1269074}{10.2307/1269074}.
#'
#' @examples
#' \dontrun{
#' library("vlad")
#' library("spcadjust")
#' set.seed(1234)
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
#' # steady state
#' RNGkind("L'Ecuyer-CMRG")
#' set.seed(12345); parallel::mc.reset.stream()
#' m <- 10^3
#' tau <- 50
#' kopt <- optimal.k(QA=2, parsonnetscores=df1$Parsonnet, coeff=coeff1)
#' eocusum.arloc.h.sim(L0=370, df=df1, k=kopt, m=m, side=1, coeff=coeff1, coeff2=coeff2, nc=nc)
#' res <- sapply(1:tau, function(i){
#'   RLS <- do.call(c, parallel::mclapply( 1:m, eocusum.adoc2.sim, k=kopt, QS=2, h=2.792265, df=df1,
#'                                         m=i, coeff=coeff1, coeff2=coeff2, side=1, mc.cores=nc))
#'   list(data.frame(cbind("ARL"=mean(RLS), "ARLSE"=sd(RLS)/sqrt(m))))
#' } )
#' a <- do.call(rbind, res)
#' plot(a[, 1], type="l")
#' }
