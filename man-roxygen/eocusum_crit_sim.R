#' @references Barnard GA (1959). Control charts and stochastic processes.
#' \emph{J R Stat Soc Series B Stat Methodol}, \strong{21}(2), pp. 239--271.
#'
#' Kemp KW (1961). The Average Run Length of the Cumulative Sum Chart
#' when a V-mask is used. \emph{J R Stat Soc Series B Stat Methodol}, \strong{23}(1),pp. 149--153.
#'
#' Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455--2473.
#'
#' @examples
#' \dontrun{
#' library(vlad)
#' library(dplyr)
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## preprocess data to 30 day mortality
#' SALL <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II")))
#' SI <- subset(SALL, phase == "I")
#' y <- subset(SALL, select = y)
#' GLM <- glm(y ~ s, data = SI, family = "binomial")
#' pi1 <- predict(GLM, type = "response", newdata = data.frame(s = SALL$s))
#' pmix <- data.frame(y, pi1, pi1)
#'
#' ## (Deterioration)
#' kopt <- optimal_k(pmix = pmix, RA = 2)
#' h <- eocusum_crit_sim(L0=370, pmix=pmix, k=kopt, side = "low", verbose=TRUE, nc=4)
#'
#' ## parameters to set up a tabular CUSUM or V-Mask (upper arm)
#' d <- h/kopt
#' theta <- atan(kopt)*180/pi
#' cbind(kopt, h, theta, d)
#'
#' ## (Improvement)
#' kopt <- optimal_k(pmix = pmix, RA = 1/2)
#' h <- eocusum_crit_sim(L0=370, pmix=pmix, k=kopt, side = "up", verbose=TRUE, nc=4)
#'
#' ## parameters to set up a tabular CUSUM or V-Mask (lower arm)
#' d <- h/kopt
#' theta <- atan(kopt)*180/pi
#' cbind(kopt, h, theta, d)
#' }
