#' @references Wittenberg P, Gan FF, Knoth S (2018).
#' “A simple signaling rule for variable life‐adjusted display derived from
#' an equivalent risk‐adjusted CUSUM chart.”
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455-2473.
#' doi: \doi{10.1002/sim.7647}.
#'
#' Taylor HM (1968). “The Economic Design of Cumulative Sum Control Charts.”
#' \emph{Technometrics}, \strong{10}(3), pp. 479-488.
#' doi: \doi{10.1080/00401706.1968.10490595}.
#'
#' Crosier R (1986). “A new two-sided cumulative quality control scheme.”
#' \emph{Technometrics}, \strong{28}(3), pp. 187-194.
#' doi: \doi{10.2307/1269074}.
#'
#' @examples
#' \dontrun{
#' data("cardiacsurgery", package = "spcadjust")
#' library("dplyr")
#'
#' ## preprocess data to 30 day mortality and subset phase I/II
#' cardiacsurgery <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II")))
#'
#' s5000 <- sample_n(cardiacsurgery, size = 5000, replace = TRUE)
#' df1 <- select(cardiacsurgery, s, y)
#' df2 <- select(s5000, s, y)
#'
#' ## estimate coefficients from logit model
#' coeff1 <- round(coef(glm(y ~ s, data = df1, family = "binomial")), 3)
#' coeff2 <- round(coef(glm(y ~ s, data = df2, family = "binomial")), 3)
#'
#' ## Number of simulation runs
#' m <- 10^3
#' ## Number of cores
#' nc <- parallel::detectCores()
#' # steady state
#' RNGkind("L'Ecuyer-CMRG")
#' m <- 10^3
#' tau <- 50
#' kopt <- optimal_k(QA = 2, df = S2I, coeff = coeff1, yemp = FALSE)
#' # eocusum_arloc_h_sim(L0 = 370, df = df1, k = kopt, m = m, side = "low", coeff = coeff1,
#'  coeff2 = coeff2, nc = nc)
#' res <- sapply(0:(tau-1), function(i){
#'   RLS <- do.call(c, parallel::mclapply( 1:m, eocusum_ad_sim, k = kopt, QS = 2, h = 2.637854,
#'   df = df1, m = i, coeff = coeff1, coeff2 = coeff2, side = "low", mc.cores = nc))
#'   list(data.frame(cbind(ARL = mean(RLS), ARLSE = sd(RLS)/sqrt(m))))
#' } )
#' RES <- data.frame(cbind(M = 0:(tau-1), do.call(rbind, res)))
#' ggplot2::qplot(x = M, y = ARL, data = RES, geom = c("line", "point")) +
#' ggplot2::theme_classic()
#' }
