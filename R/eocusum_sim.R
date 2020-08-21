#' @name optimal_k
#' @title Compute approximate optimal k
#' @description Compute approximate optimal k.
#'
#' @param pmix Data Frame. A three column data frame. First column is the operation outcome.
#' Second column are the predicted probabilities from the risk model. Third column can be either the
#'  predicted probabilities from the risk model or average outcome.
#' @param RA Double. Odds ratio of death under the alternative hypotheses. Detecting deterioration
#' in performance with increased mortality risk by doubling the odds Ratio \code{RA = 2}. Detecting
#'  improvement in performance with decreased mortality risk by halving the odds ratio of death
#'  \code{RA = 1/2}. Odds ratio of death under the null hypotheses is \code{1}.
#' @param yemp Logical. If \code{TRUE}, use emirical outcome values, else use model.
#'
#' @return Returns a single value which is the approximate optimal \code{k}.
#'
#' @details Formula deterioration:  \deqn{ k{det} = \frac{R{A} - 1 - log(R{A})}{log(R{A})}\bar{p} ,
#' R{A} > 1    }
#'          Formula improvement:    \deqn{ k{imp} = \frac{1 - R{A} + log(R{A})}{log(R{A})}\bar{p} ,
#' R{A} < 1    }
#'
#' @template optimal_k
#'
#' @author Philipp Wittenberg
#' @export
optimal_k <- function(pmix, RA, yemp = FALSE) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(RA, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_logical(yemp, len = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .optimal_k(pmix, RA, yemp)
}

#' @name eocusum_arl_sim
#' @title Compute ARLs of EO-CUSUM control charts using simulation
#' @description Compute ARLs of EO-CUSUM control charts using simulation.
#'
#' @param pmix Data Frame. A three column data frame. First column is the operation outcome.
#' Second column are the predicted probabilities from the risk model. Third column can be either the
#'  predicted probabilities from the risk model or average outcome.
#' @param r Integer. Number of of simulation runs.
#' @param k Double. Reference value of the CUSUM control chart. Either \code{0} or a positive
#' value. Can be determined with function \code{\link{optimal_k}}.
#' @param h Double. Decision interval (alarm limit, threshold) of the CUSUM control chart.
#' @param RQ Double. Defines the true performance of a surgeon with the odds ratio ratio of death
#' \code{RQ}. Use \code{RQ = 1} to compute the in-control ARL and other values to compute the
#' out-of-control ARL.
#' @param yemp Logical. If \code{TRUE} use observed outcome value, if \code{FALSE} use estimated
#' binary logistc regression model.
#' @param side Character. Default is \code{"low"} to calculate ARL for the upper arm of the V-mask.
#'  If side = \code{"up"}, calculate the lower arm of the V-mask.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template eocusum_arl_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_arl_sim <- function(r, pmix, k, h, RQ = 1, yemp = FALSE, side = "low") {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_integerish(r, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(k, len = 1, add = arg_checks)
  checkmate::assert_numeric(h, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_logical(yemp, len = 1, add = arg_checks)
  side <- tolower(side)
  checkmate::assert_choice(side, choices = c("low", "up"), add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  iside <- switch(side, low = 1, up = 2)
  if (RQ < 1 && iside == 1) stop("For detecting deterioration (side='low') RQ must a positive numeric value >= 1")
  if (RQ > 1 && iside == 2) stop("For detecting improvement (side='up') RQ must a positive numeric value <= 1")
  .eocusum_arl_sim(r, pmix, k, h, RQ, yemp, iside)
}

#' @name eocusum_ad_sim
#' @title Compute steady-state ARLs of EO-CUSUM control charts using simulation
#' @description Compute steady-state ARLs of EO-CUSUM control charts using simulation.
#'
#' @inheritParams eocusum_arl_sim
#' @param m Integer. Simulated in-control observations.
#' @param type Character. Default argument is \code{"cond"} for computation of conditional
#' steady-state. Other option is the cyclical steady-state \code{"cycl"}.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template eocusum_ad_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_ad_sim <- function(r, pmix, k, h, RQ = 1, side = "low", type = "cond", m = 50) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_integerish(r, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(k, len = 1, add = arg_checks)
  checkmate::assert_numeric(h, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  side <- tolower(side)
  checkmate::assert_choice(side, choices = c("low", "up"), add = arg_checks)
  iside <- switch(as.character(side), low = 1, up = 2)
  if (RQ < 1 && iside == 1) stop("For detecting deterioration (side='low') RQ must a positive numeric value >= 1")
  if (RQ > 1 && iside == 2) stop("For detecting improvement (side='up') RQ must a positive numeric value <= 1")
  type <- tolower(type)
  checkmate::assert_choice(type, choices = c("cond", "cycl"), add = arg_checks)
  itype <- switch(type, cond = 1, cycl = 2)
  checkmate::assert_integerish(m, lower = 1, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  .eocusum_ad_sim(r, pmix, k, h, RQ, iside, itype, m)
}

#' @name eocusum_crit_sim
#' @title Compute alarm threshold of EO-CUSUM control charts using simulation
#' @description Compute alarm threshold of EO-CUSUM control charts using simulation.
#'
#' @param L0 Double. Prespecified in-control Average Run Length.
#' @inheritParams eocusum_arl_sim
#' @param m Integer. Number of simulation runs.
#' @param nc Integer. Number of cores used for parallel processing. Value is passed to
#' \code{\link{parSapply}}.
#' @param hmax Integer. Maximum value of \code{h} for the grid search.
#' @param jmax Integer. Number of digits for grid search.
#' @param verbose Logical. If \code{TRUE} verbose output is included, if \code{FALSE} a quiet
#' calculation of \code{h} is done.
#'
#' @return Returns a single value which is the control limit \code{h} for a given in-control ARL.
#'
#' @template eocusum_crit_sim
#'
#' @details Determines the control limit ("\code{h}") for given in-control ARL (\code{"L0"})
#' applying a grid search using \code{\link{eocusum_arl_sim}} and \code{\link{parSapply}}.
#'
#' @author Philipp Wittenberg
#'
#' @export
eocusum_crit_sim <- function(L0, pmix, k, RQ = 1, side = "low", yemp = FALSE, m = 1e4, nc = 1, hmax = 30, jmax = 4, verbose = FALSE) {
  arg_checks <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(L0, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_data_frame(pmix, ncols = 3, add = arg_checks)
  checkmate::assert_numeric(k, len = 1, add = arg_checks)
  checkmate::assert_numeric(RQ, len = 1, lower = 0, add = arg_checks)
  checkmate::assert_logical(yemp, len = 1, add = arg_checks)
  checkmate::assert_integerish(m, lower = 1, add = arg_checks)
  checkmate::assert_integerish(nc, lower = 1, add = arg_checks)
  checkmate::assert_integerish(hmax, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_integerish(jmax, len = 1, lower = 1, add = arg_checks)
  checkmate::assert_logical(verbose, len = 1, add = arg_checks)
  side <- tolower(side)
  checkmate::assert_choice(side, choices = c("low", "up"), add = arg_checks)
  # iside <- switch(side, low = 1, up = 2)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  # iside <- switch(as.character(side), low = 1, up = 2)
  if (RQ < 1 && side == "low") stop("For detecting deterioration (side='low') RQ must a positive numeric value >= 1")
  if (RQ > 1 && side == "up") stop("For detecting improvement (side='up') RQ must a positive numeric value <= 1")

  cl <- parallel::makeCluster(getOption("cl.cores", nc))

  for ( h in 1:hmax ) {
    parallel::clusterExport(cl, c("h", "eocusum_arl_sim", "pmix", "k", "RQ", "side", "yemp"), envir = environment())
    # iside <- switch(iside, 1 = "low", 2 = "up")
    L1 <- mean(parallel::parSapply(cl, 1:m, eocusum_arl_sim, pmix = pmix, k = k, h = h, RQ = RQ, side = side, yemp = yemp))
    if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
    if ( L1 > L0 ) break
  }
  h1 <- h

  for ( j in 1:jmax ) {
    for ( dh in 1:19 ) {
      h <- h1 + (-1)^j*dh/10^j
      parallel::clusterExport(cl, c("h", "eocusum_arl_sim", "pmix", "k", "RQ", "side", "yemp"), envir = environment())
      # iside <- switch(iside, 1 = "low", 2 = "up")
      L1 <- mean(parallel::parSapply(cl, 1:m, eocusum_arl_sim, pmix = pmix, k = k, h = h, RQ = RQ, side = side, yemp = yemp))
      if ( verbose ) cat(paste("h =", h, "\tARL =", L1, "\n"))
      if ( (j %% 2 == 1 & L1 < L0) | (j %% 2 == 0 & L1 > L0) ) break
    }
    h1 <- h
  }

  if ( L1 < L0 ) h <- h + 1/10^jmax
  parallel::stopCluster(cl)
  h
}
