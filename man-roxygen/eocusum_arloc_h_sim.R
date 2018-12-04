#' @references Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455-2473.
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
#'
#' ## Lower CUSUM (detecting deterioration)
#' ## k = 0
#' eocusum_arloc_h_sim(L0 = 370, df = df1, k = 0, m = m, side = "low", coeff = coeff1, coeff2 =
#' coeff2, nc = nc)
#' ## use function optimal_k() to determine k = kopt
#' kopt <- optimal_k(QA = 2, df = S2I, coeff = coeff1, yemp = FALSE)
#' eocusum_arloc_h_sim(L0 = 370, df = df1, k = kopt, m = m, side = "low", coeff = coeff1, coeff2 =
#' coeff2, nc = nc)
#'
#' ## Upper CUSUM (detecting improvement)
#' ## k = 0
#' eocusum_arloc_h_sim(L0 = 370, df = df1, k = 0, m = m, side = "up", coeff = coeff1, coeff2 =
#' coeff2, nc = nc)
#' ## use function optimal_k() to determine k = kopt
#' kopt <- optimal_k(QA = 1/2, df = S2I, coeff = coeff1, yemp = FALSE)
#' eocusum_arloc_h_sim(L0 = 370, df = df1, k = kopt, m = m, side = "up", coeff = coeff1, coeff2 =
#'  coeff2, nc = nc)
#' }
