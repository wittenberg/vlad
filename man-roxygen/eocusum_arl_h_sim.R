#' @references Barnard GA (1959). "Control charts and stochastic processes."
#' \emph{J R Stat Soc Series B Stat Methodol}, \strong{21}(2), pp. 239-271.
#'
#' Kemp KW (1961). "The Average Run Length of the Cumulative Sum Chart
#' when a V-mask is used." \emph{J R Stat Soc Series B Stat Methodol}, \strong{23}(1),pp. 149-153.
#' doi: \href{https://doi.org/10.2307/2985287}{10.2307/2985287}.
#'
#' @examples
#' \dontrun{
#' library("vlad")
#' library("spcadjust")
#' set.seed(1234)
#' data("cardiacsurgery")
#' df1 <- subset(cardiacsurgery, select=c(Parsonnet, status))
#' ## estimate coefficients from logit model
#' coeff1 <- round(coef(glm(status~Parsonnet, data=df1, family="binomial")), 3)
#' ## Number of simulation runs
#' m <- 10^3
#' ## Number of cores
#' nc <- parallel::detectCores()
#' ## Detect deterioration
#' QA <- 2
#' kopt <- optimal_k(QA=QA, df=df1, coeff=coeff1)
#' h <- eocusum_arl_h_sim(L0=370, df=df1, k=kopt, m=m, coeff=coeff1, side="low", nc=nc)
#' ## V-Mask parameters
#' d <- h/kopt
#' theta <- atan(kopt)*180/pi
#' cbind(kopt, h, theta, d)
#' }
