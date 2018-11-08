#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  “Monitoring surgical performance using risk-adjusted cumulative sum charts.”
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441-452.
#'  doi: \doi{10.1093/biostatistics/1.4.441}.
#'
#' Brook D and Evans DA (1972)
#'  "An approach to the probability distribution of CUSUM run length"
#'  \emph{Biometrika}, \strong{59}(3), pp. 539-549
#'  doi: \doi{10.1093/biomet/59.3.539}.
#'
#' Webster RA and Pettitt AN (2007)
#'  \emph{Communications in Statistics - Simulation and Computation} \strong{36}(3), pp. 471-482
#'  doi: \doi{10.1080/03610910701208361}.
#'
#' @examples
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
