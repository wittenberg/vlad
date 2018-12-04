#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  “Monitoring surgical performance using risk-adjusted cumulative sum charts.”
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441-452.
#'
#' Wittenberg P, Gan FF, Knoth S (2018).
#' “A simple signaling rule for variable life‐adjusted display derived from
#' an equivalent risk‐adjusted CUSUM chart.”
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455-2473.
#'
#' @examples
#' \dontrun{
#' library("dplyr")
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## preprocess data to 30 day mortality and subset phase I (In-control) of surgeons 2
#' S2I <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
#'   filter(phase == "I", surgeon == 2) %>% select(s, y)
#'
#' ## estimate coefficients from logit model
#' coeff1 <- round(coef(glm(y ~ s, data = S2I, family = "binomial")), 3)
#'
#' ## control limit for detecting deterioration RA = 2:
#' racusum_crit_sim(L0 = 740, df = S2I, coeff = coeff1, m = 10^3, nc = 4)
#'}
