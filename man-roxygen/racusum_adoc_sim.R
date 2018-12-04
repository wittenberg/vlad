#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#' Monitoring surgical performance using risk-adjusted cumulative sum charts.
#' \emph{Biostatistics}, \strong{1}(4), pp. 441-452.
#' doi: \doi{10.1093/biostatistics/1.4.441}.
#'
#' Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455-2473.
#'
#' Taylor HM (1968). The Economic Design of Cumulative Sum Control Charts.
#' \emph{Technometrics}, \strong{10}(3), pp. 479-488.
#'
#' Crosier R (1986). A new two-sided cumulative quality control scheme.
#' \emph{Technometrics}, \strong{28}(3), pp. 187-194.
#'
#' @examples
#' \dontrun{
#' data("cardiacsurgery", package="spcadjust")
#' # build data set
#' df1 <- subset(cardiacsurgery, select=c(Parsonnet, status))
#'
#' # estimate coefficients from logit model
#' coeff1 <- round(coef(glm(status ~ Parsonnet, data=df1, family="binomial")), 3)
#'
#' # simulation of conditional steady state
#' m <- 10^3
#' tau <- 50
#' res <- sapply(0:(tau-1), function(i){
#'  RLS <- do.call(c, parallel::mclapply( 1:m, racusum_adoc_sim, RQ=2, h=2.0353, df=df1, m=i,
#'                                        coeff=coeff1, coeff2=coeff1,
#'                                        mc.cores=parallel::detectCores()) )
#'  list(data.frame(cbind(ARL=mean(RLS), ARLSE=sd(RLS)/sqrt(m))))
#' } )
#'
#' # plot
#' RES <- data.frame(cbind(M=0:(tau-1), do.call(rbind, res)))
#' ggplot2::qplot(x=M, y=ARL, data=RES, geom=c("line", "point")) +
#' ggplot2::theme_classic()
#' }
