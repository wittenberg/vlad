#' @references Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455--2473.
#'
#' @examples
#' library("dplyr")
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## preprocess data to 30 day mortality and subset phase I (In-control) of surgeons 2
#' S2I <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'         phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
#'   filter(phase == "I", surgeon == 2) %>% select(s, y)
#'
#' coeff <- coef(glm(y ~ s, data = S2I, family = "binomial"))
#'
#' ## (Deterioration)
#' optimal_k(QA = 2, df = S2I, coeff = coeff, yemp = FALSE)
#'
#' ## manually find optimal k for detecting deterioration
#' QA <- 2
#' pbar <- mean(sapply(S2I[, 1], gettherisk, coef = coeff))
#' kopt <- pbar * ( QA - 1 - log(QA) ) / log(QA)
#'
#' all.equal(kopt, optimal_k(QA = 2, df = S2I, coeff = coeff, yemp = FALSE))
#'
#' ## (Improvement)
#' optimal_k(QA = 1/2, df = S2I, coeff = coeff, yemp = FALSE)
#'
#' ## manually find optimal k for detecting improvement
#' QA <- 1/2
#' pbar <- mean(sapply(S2I[, 1], gettherisk, coef = coeff))
#' kopt <- pbar * ( 1 - QA + log(QA) ) / log(QA)
#'
#' all.equal(kopt, optimal_k(QA = 1/2, df = S2I, coeff = coeff, yemp = FALSE))
