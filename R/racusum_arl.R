#' @name racusum_arl
#' @title ARL of RA-CUSUM charts
#' @description Compute the ARL of risk-adjusted CUSUM charts.
#' @param pmix Numeric Matrix. A three column matrix. First column is the risk
#' score distribution. Second column are the predicted probabilities from the risk model. Third
#'  column can be either the predicted probabilities from the risk model or average outcome per
#'  risk score, see examples.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}. Odds ratio of death under the null hypotheses is \code{1}.
#' \code{RQ}. Use \code{RQ = 1} to compute the in-control ARL and other values to compute the
#' out-of-control ARL.
#' @param RQ Double. Defines the true performance of a surgeon with the odds ratio ratio of death
#' \code{RQ}. Use \code{RQ = 1} to compute the in-control ARL and other values to compute the
#' out-of-control ARL.
#' @param h Double. \code{h} is the control limit (>\code{0}).
#' @param scaling Double. The \code{scaling} parameter controls the quality of the approximation,
#' larger values achieve higher accuracy but increase the computation burden (larger transition
#' probability matrix).
#' @param rounding Character. If \code{rounding = "p"} a paired rounding implementation of
#' \emph{Knoth et al. (2019)} is used, if \code{rounding = "s"} a simple rounding method of
#' \emph{Steiner et al. (2000)} is used.
#' @param method Character. If \code{method = "Toep"} a combination of Sequential Probability Ratio
#'  Test and Toeplitz matrix structure is used to calculate the ARL. \code{"ToepInv"} computes the
#'  inverted matrix using Toeplitz matrix structure. \code{"BE"} solves a linear equation system
#'  using the classical approach of \emph{Brook and Evans (1972)} to calculate the ARL.
#' @return Returns a single value which is the Average Run Length for \code{"racusum_arl_mc"} and
#' the Run Length for \code{"racusum_arl_sim"}.
#' @author Philipp Wittenberg
#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#'  Monitoring surgical performance using risk-adjusted cumulative sum charts.
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441--452.
#'
#' Knoth S, Wittenberg P and Gan FF (2019).
#' Risk-adjusted CUSUM charts under model error.
#' \emph{Statistics in Medicine}, \strong{38}(12), pp. 2206--2218.
#'
#' Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455--2473.
#'
#' Brook D and Evans DA (1972)
#'  An approach to the probability distribution of CUSUM run length.
#'  \emph{Biometrika}, \strong{59}(3), pp. 539--549
#'
#' Webster RA and Pettitt AN (2007)
#' Stability of approximations of average run length of risk-adjusted CUSUM schemes using
#' the Markov approach: comparing two methods of calculating transition probabilities.
#'  \emph{Communications in Statistics - Simulation and Computation} \strong{36}(3), pp. 471--482
#' @examples
#' \dontrun{
#' library(vlad)
#' library(dplyr)
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## Markov Chain
#' ## preprocess data to 30 day mortality and subset phase I (In-control) of surgeons 2
#' SALLI <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'         phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
#'   filter(phase == "I") %>% select(s, y)
#'
#' ## estimate risk model, get relative frequences and probabilities
#' mod1 <- glm(y ~ s, data = SALLI, family = "binomial")
#' fi  <- as.numeric(table(SALLI$s) / length(SALLI$s))
#' usi <- sort(unique(SALLI$s))
#' pi1 <- predict(mod1, newdata = data.frame(s = usi), type = "response")
#' pi2 <- tapply(SALLI$y, SALLI$s, mean)
#'
#' ## set up patient mix (risk model)
#' pmix1  <- data.frame(fi, pi1, pi1)
#'
#' ## Average Run Length for detecting deterioration RA = 2:
#' racusum_arl_mc(pmix = pmix1, RA = 2, RQ = 1, h = 4.5)
#'
#' ## Average Run Length for detecting improvement RA = 1/2:
#' racusum_arl_mc(pmix = pmix1, RA = 1/2, RQ = 1, h = 4)
#'
#' ## set up patient mix (model free)
#' pmix2  <- data.frame(fi, pi1, pi2)
#'
#' ## Average Run Length for detecting deterioration RA = 2:
#' racusum_arl_mc(pmix = pmix2, RA = 2, RQ = 1, h = 4.5)
#'
#' ## Average Run Length for detecting improvement RA = 1/2:
#' racusum_arl_mc(pmix = pmix2, RA = 1/2, RQ = 1, h = 4)
#'
#' ## compare results with R-code function 'findarl()' from Steiner et al. (2000)
#' source("https://bit.ly/2KC0SYD")
#' all.equal(findarl(pmix = pmix1, R1 = 2, R = 1, CL = 4.5, scaling = 600),
#'          racusum_arl_mc(pmix = pmix1, RA = 2, RQ = 1, h = 4.5, scaling = 600, rounding = "s"))
#'
#' ## Monte Carlo simulation
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
#' @export
racusum_arl_mc <- function(h, pmix, RA, RQ, scaling = 600, rounding = "p", method = "Toep") {
  arg_checks <- checkmate::makeAssertCollection()
  pmix <- as.matrix(pmix)
  checkmate::assert_matrix(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(scaling, lower = 1, add = arg_checks)
  rounding <- tolower(rounding)
  checkmate::assert_choice(rounding, c("p", "s"))
  irounding <- switch(rounding, p = 1, s = 2)
  checkmate::assert_choice(method, c("Toep", "ToepInv", "BE"))
  imethod <- switch(method, Toep = 1, ToepInv = 2, BE = 3)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_arl_mc(pmix, RA, RQ, h, scaling, irounding, imethod)
}

#' @rdname racusum_arl
#' @param r Integer. Number of runs.
#' @param pmix Data Frame. A three column data frame. First column is the operation outcome.
#' Second column are the predicted probabilities from the risk model. Third column can be either the
#'  predicted probabilities from the risk model or average outcome.
#' @param yemp Logical. If \code{TRUE} use observed outcome value, if \code{FALSE} use estimated
#' binary logistc regression model.
#' @export
racusum_arl_sim <- function(h, pmix, r, RA = 2, RQ = 1, yemp = FALSE) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_integerish(r, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(h, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_logical(yemp, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_arl_sim(r, pmix, h, RA, RQ, yemp)
}
