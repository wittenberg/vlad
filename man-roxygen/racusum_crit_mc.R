#' @references Knoth S, Wittenberg P and Gan FF (2019).
#' Risk-adjusted CUSUM charts under model error.
#' \emph{Statistics in Medicine}, \strong{38}(12), pp. 2206--2218.
#'
#' Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  Monitoring surgical performance using risk-adjusted cumulative sum charts.
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441--452.
#'
#' Brook D and Evans DA (1972)
#'  An approach to the probability distribution of CUSUM run length.
#'  \emph{Biometrika}, \strong{59}(3), pp. 539--549
#'
#' Webster RA and Pettitt AN (2007)
#' Stability of approximations of average run length of risk-adjusted CUSUM schemes using
#' the Markov approach: comparing two methods of calculating transition probabilities.
#'  \emph{Communications in Statistics - Simulation and Computation} \strong{36}(3), pp. 471--482
#'
#' @examples
#' \dontrun{
#' library(vlad)
#' library(dplyr)
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## preprocess data to 30 day mortality and subset phase I (In-control) of surgeons 2
#' S2I <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
#'   filter(phase == "I", surgeon == 2) %>% select(s, y)
#'
#' ## estimate risk model, get relative frequences and probabilities
#' mod1 <- glm(y ~ s, data = S2I, family = "binomial")
#' fi  <- as.numeric(table(S2I$s) / length(S2I$s))
#' usi <- sort(unique(S2I$s))
#' pi1 <- predict(mod1, newdata = data.frame(s = usi), type = "response")
#'
#' ## set up patient mix
#' pmix  <- data.frame(fi, pi1, pi1)
#'
#' ## control limit for detecting deterioration RA = 2:
#' racusum_crit_mc(pmix = pmix, L0 = 740, RA = 2, RQ = 1)
#' ## control limit for detecting improvement RA = 1/2:
#' racusum_crit_mc(pmix = pmix, L0 = 740, RA = 0.5, RQ = 1)
#' }
