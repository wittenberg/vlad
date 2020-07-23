#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  Monitoring surgical performance using risk-adjusted cumulative sum charts.
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441--452.
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
#' set.seed(1234)
#' SALLI <- cardiacsurgery %>% mutate(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
#'   filter(phase == "I") %>% select(s, y)
#'
#' ## estimate risk model, get relative frequences and probabilities
#' mod1 <- glm(y ~ s, data = SALLI, family = "binomial")
#' y <- SALLI$y
#' pi1 <- fitted.values(mod1)
#'
#' ## set up patient mix (risk model)
#' pmix <- data.frame(y, pi1, pi1)
#' h <- 2.75599
#'
#' m <- 1e4
#' RLS <- sapply(1:m, racusum_arl_sim, h=h, pmix=pmix, RA=2)
#' data.frame(cbind(ARL=mean(RLS), ARLSE=sd(RLS)/sqrt(m), h, m))
#' }
