#' @references Barnard GA (1959). "Control charts and stochastic processes."
#' \emph{J R Stat Soc Series B Stat Methodol}, \strong{21}(2), pp. 239-271.
#'
#' Kemp KW (1961). "The Average Run Length of the Cumulative Sum Chart
#' when a V-mask is used." \emph{J R Stat Soc Series B Stat Methodol}, \strong{23}(1),pp. 149-153.
#' doi: \href{https://doi.org/10.2307/2985287}{10.2307/2985287}.
#'
#' Wittenberg P, Gan FF, Knoth S (2018).
#' “A simple signaling rule for variable life‐adjusted display derived from
#' an equivalent risk‐adjusted CUSUM chart.”
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455-2473.
#' doi: \href{https://doi.org/10.1002/sim.7647}{10.1002/sim.7647}.
#'
#' @examples
#' \dontrun{
#' data("cardiacsurgery", package = "spcadjust")
#' library("dplyr")
#'
#' ## preprocess data to 30 day mortality and subset phase I (In-control) of surgeons 2
#' S2I <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
#'   filter(phase == "I", surgeon == 2) %>% select(s, y)
#'
#' ## estimate coefficients from logit model
#' coeff1 <- coef(glm(y ~ s, data = S2I, family = "binomial"))
#' ## Number of simulation runs
#' m <- 10^3
#' set.seed(1234)
#' ## Number of cores
#' nc <- parallel::detectCores()
#'
#' ## determine k for detecting deterioration
#' kopt <- optimal_k(QA = 2, df = S2I, coeff = coeff, yemp = FALSE)
#'
#' ## compute threshold for prespecified in-control ARL
#' h <- eocusum_crit_sim(L0 = 370, df = S2I, k = kopt, m = m, coeff = coeff1, side = "low",
#' nc = nc)
#'
#' ## parameters to set up a tabular CUSUM or V-Mask
#' d <- h/kopt
#' theta <- atan(kopt)*180/pi
#' cbind(kopt, h, theta, d)
#' }
