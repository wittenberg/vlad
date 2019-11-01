#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  Monitoring surgical performance using risk-adjusted cumulative sum charts.
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441--452.
#'  \doi{10.1093/biostatistics/1.4.441}
#'
#' Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455--2473.
#' \doi{10.1002/sim.7647}
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
#' h <- racusum_crit_sim(pmix = pmix, L0 = 370, RA = 2, nc = 4, verbose = TRUE)
#' }
