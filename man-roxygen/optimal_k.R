#' @references Wittenberg P, Gan FF, Knoth S (2018).
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
#' optimal_k(pmix = pmix, RA = 2, yemp = FALSE)
#'
#' ## manually find optimal k for detecting deterioration
#' RA <- 2
#' pbar <- mean(pmix$pi1)
#' kopt <- pbar * ( RA - 1 - log(RA) ) / log(RA)
#'
#' all.equal(kopt, optimal_k(pmix = pmix, RA = 2, yemp = FALSE))
#'
#' ## (Improvement)
#' optimal_k(pmix = pmix, RA = 1/2, yemp = FALSE)
#'
#' ## manually find optimal k for detecting improvement
#' RA <- 1/2
#' pbar <- mean(pmix$pi1)
#' kopt <- pbar * ( 1 - RA + log(RA) ) / log(RA)
#'
#' all.equal(kopt, optimal_k(pmix = pmix, RA = 1/2, yemp = FALSE))
#' }
