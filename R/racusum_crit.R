#' @name racusum_crit
#' @title Alarm thresholds of RA-CUSUM charts
#' @description Compute alarm threshold of risk-adjusted CUSUM charts.
#' @param L0 Double. Prespecified Average Run Length.
#' @param jmax Integer. Number of digits for grid search.
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
#' @param verbose Logical. If \code{FALSE} a quiet calculation of \code{h} is done. If \code{TRUE}
#'  verbose output of the search procedure is included.
#' @return Returns a single value which is the control limit \code{h} for a given In-control ARL.
#' @details
#' Determines the control limit for given in-control ARL (\code{"L0"}) using
#' \code{\link{racusum_arl_mc}} by applying a grid search.
#' @author Philipp Wittenberg
#' @references Knoth S, Wittenberg P and Gan FF (2019).
#' Risk-adjusted CUSUM charts under model error.
#' \emph{Statistics in Medicine}, \strong{38}(12), pp. 2206--2218.
#'
#' Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455--2473.
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
#' @examples
#' \dontrun{
#' library(vlad)
#' library(dplyr)
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## Markov Chain
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
#'
#' ## Monte Carlo simulation
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
#' @export
racusum_crit_mc <- function(L0, pmix, RA, RQ, scaling = 600, rounding = "p", method = "Toep", jmax = 4, verbose = FALSE) {
  arg_checks <- checkmate::makeAssertCollection()
  pmix <- as.matrix(pmix)
  checkmate::assert_matrix(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_integerish(L0, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_integerish(scaling, lower = 1, add = arg_checks)
  rounding <- tolower(rounding)
  checkmate::assert_choice(rounding, c("p", "s"))
  irounding <- switch(rounding, p = 1, s = 2)
  checkmate::assert_choice(method, c("Toep", "ToepInv", "BE"))
  imethod <- switch(method, Toep = 1, ToepInv = 2, BE = 3)
  checkmate::assert_integerish(jmax, lower = 0, add = arg_checks)
  checkmate::assert_logical(verbose, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .racusum_crit_mc(pmix, L0, RA, RQ, scaling, irounding, imethod, jmax, verbose)
}

#' @rdname racusum_crit
#' @param yemp Logical. If \code{TRUE}, use emirical outcome values, else use model.
#' @param m Integer. Number of simulation runs.
#' @param nc Integer. Number of cores used for parallel processing. Value is passed to
#' \code{\link{parSapply}}.
#' @param hmax Integer. Maximum value of \code{h} for the grid search.
#' @details Determines the control limit ("\code{h}") for given in-control ARL (\code{"L0"})
#' applying a grid search using \code{\link{racusum_arl_sim}} and \code{\link{parSapply}}.
#' @export
racusum_crit_sim <- function(L0, pmix, RA = 2, RQ = 1, yemp = FALSE, m = 1e4, nc = 1, hmax = 30, jmax = 4, verbose = FALSE) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(L0, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_logical(yemp, len = 1, add = arg_checks)
  checkmate::assert_integerish(m, lower = 1, add = arg_checks)
  checkmate::assert_integerish(nc, lower = 1, add = arg_checks)
  checkmate::assert_integerish(hmax, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_integerish(jmax, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_logical(verbose, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)

  cl <- parallel::makeCluster(getOption("cl.cores", nc))

  for ( h in 1:hmax ) {
    parallel::clusterExport(cl, c("h", "racusum_arl_sim", "pmix", "RA", "RQ", "yemp", "m"), envir = environment())
    L1 <- mean(parallel::parSapply(cl, 1:m, racusum_arl_sim, pmix = pmix, h = h, RA = RA, RQ = RQ, yemp = yemp))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      parallel::clusterExport(cl, c("h", "racusum_arl_sim", "pmix", "RA", "RQ", "yemp", "m"), envir = environment())
      L1 <- mean(parallel::parSapply(cl, 1:m, racusum_arl_sim, pmix = pmix, h = h, RA = RA, RQ = RQ, yemp = yemp))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax
  parallel::stopCluster(cl)
  h
}

