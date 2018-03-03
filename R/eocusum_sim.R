#' @name optimal_k
#' @title Compute optimal k
#' @description Compute optimal k.
#'
#' @param QA double. Defines the performance of a surgeon with the odds ratio ratio of death Q.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model. For more information see details.
#' @param parsonnetscores NumericVector. Vector of Parsonnet Scores.
#'
#' @return Returns a single value which is the approximate optimal k for a set of given Parsonnet scores.
#'
#' @details Formula deterioration:  \deqn{ k{det} = \frac{Q{A} - 1 - log(Q{A})}{log(Q{A})}\bar{p} , Q{A} > 1    }
#'          Formula improvement:    \deqn{ k{imp} = \frac{1 - Q{A} + log(Q{A})}{log(Q{A})}\bar{p} , Q{A} < 1    }
#'
#' @template optimal_k
#'
#' @author Philipp Wittenberg
#' @export
optimal_k <- function(QA, parsonnetscores, coeff) {
  if (is.null(QA) || is.na(QA) || QA <= 0) {stop("QA must a positive numeric value")}
  QA <- as.numeric(QA)
  if (is.null(parsonnetscores) || is.na(parsonnetscores) || is.vector(parsonnetscores) != "TRUE") {stop("Argument 'parsonnetscore' must be an integer value")}
  parsonnetscores <- as.vector(parsonnetscores)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  .optimal_k(QA, parsonnetscores, coeff)
}

#' @name gettherisk
#' @title Compute Risk of death
#' @description Compute Risk of death.
#'
#' @param parsonnetscore int. Parsonnet Score.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model.
#'
#' @return Returns a single value which is the expected risk based on a risk model.
#'
#' @template gettherisk
#'
#' @author Philipp Wittenberg
#'
#' @export
gettherisk <- function(parsonnetscore, coeff) {
  if (is.null(parsonnetscore) || is.na(parsonnetscore) || is.integer(parsonnetscore) != "TRUE") {stop("Argument 'parsonnetscore' must be an integer value")}
  parsonnetscore <- as.integer(parsonnetscore)
  if (is.null(coeff) || is.na(coeff) || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  .gettherisk(parsonnetscore, coeff)
}

#' @name calceo
#' @title Compute Expected minus Observed value
#' @description Compute Expected minus Observed value.
#'
#' @param df DataFrame. First column Parsonnet Score and second column outcome of each operation.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model.
#' @param yemp boolean. If TRUE use observed outcome value, if FALSE use estimated binary logistc
#'  regression model.
#'
#' @return Returns a single value which is the difference between expected risk and observed outcome.
#'
#' @template calceo
#'
#' @keywords keyword1 keyword2
#' @family VLAD functions
#' @author Philipp Wittenberg
#' @export
calceo <- function(df, coeff, yemp = TRUE) {
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff) || length(coeff) != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.na(yemp) || is.logical(yemp) != "TRUE") {warning("Argument 'yemp' must be logical using TRUE as default value")}
  yemp <- as.logical(yemp)
  .calceo(df, coeff, yemp)
}

#' @name eocusum_arl_sim
#' @title Compute ARLs of EO-CUSUM control charts using simulation
#' @description Compute ARLs of EO-CUSUM control charts using simulation.
#'
#' @param r int. Number of of simulation runs.
#' @param k double. Reference value of the CUSUM control chart.
#' @param h double. Decision interval (alarm limit, threshold) of the CUSUM control chart.
#' @inheritParams calceo
#' @param side character. Default is "low" to calculate ARL for the upper arm of the V-mask. If side = "up",
#'  calculate the lower arm of the V-mask.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template eocusum_arl_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_arl_sim <- function(r, k, h, df, coeff, yemp = TRUE, side = "low") {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff) || length(coeff) != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.na(yemp) || is.logical(yemp) != "TRUE") {warning("Argument 'yemp' must be logical using TRUE as default value")}
  yemp <- as.logical(yemp)
  iside <- switch(as.character(side), low = 1, up = 2)
  if (is.null(iside)) {
    warning("No valid input, u sing side='low' (deterioration) as default")
    iside <- 1
  }
  .eocusum_arl_sim(r, k, h, df, coeff, yemp, iside)
}

#' @name eocusum_arloc_sim
#' @title Compute Out of Control ARLs of EO-CUSUM control charts using simulation
#' @description Compute Out of Control ARLs of EO-CUSUM control charts using simulation.
#'
#' @inheritParams eocusum_arl_sim
#' @param coeff2 NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model of a resampled dataset.
#' @param QS double. Defines the performance of a surgeon with the odds ratio ratio of death Q.
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template eocusum_arloc_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_arloc_sim <- function(r, k, h, df, coeff, coeff2, QS = 1, side = "low") {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.null(coeff2) || is.na(coeff2) || length(coeff2) != 2) {stop("Model coefficients 'coeff2' must be a numeric vector with two elements")}
  coeff2 <- as.vector(coeff2)
  iside <- switch(as.character(side), low = 1, up = 2)
  if (is.null(iside)) {
    warning("No valid input, u sing side='low' (deterioration) as default")
    iside <- 1
  }
  QS <- as.numeric(QS)
  if (is.na(QS) || QS < 0) {stop("QS must a positive numeric value")}
  else if (QS < 1 && iside == 1) {stop("For detecting deterioration (side='low') QS must a positive numeric value >= 1")}
  else if (QS > 1 && iside == 2) {stop("For detecting improvement (side='up') QS must a positive numeric value <= 1")}
  .eocusum_arloc_sim(r, k, h, df, coeff, coeff2, QS, iside)
}

#' @name eocusum_adoc_sim
#' @title Compute steady-state ARLs of EO-CUSUM control charts using simulation
#' @description Compute steady-state ARLs of EO-CUSUM control charts using simulation.
#'
#' @inheritParams eocusum_arloc_sim
#' @param m Integer. Simulated in-control observations.
#' @param type character. Default argument is "cond" for computation of conditional steady-state.
#' Other option is the cyclical steady-state "cycl".
#'
#' @return Returns a single value which is the Run Length.
#'
#' @template eocusum_adoc_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_adoc_sim <- function(r, k, h, df, coeff, coeff2, QS = 1, side = "low", type = "cond", m = 50) {
  r <- as.integer(r)
  if (is.na(r) || r <= 0) {stop("Number of simulation runs 'r' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  h <- as.numeric(h)
  if (is.na(h) || h <= 0) {stop("Control limit 'h' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.null(coeff2) || is.na(coeff2) || length(coeff2) != 2) {stop("Model coefficients 'coeff2' must be a numeric vector with two elements")}
  coeff2 <- as.vector(coeff2)
  iside <- switch(as.character(side), low = 1, up = 2)
  if (is.null(iside)) {
    warning("No valid input, using side='low' (deterioration) as default")
    iside <- 1
  }
  QS <- as.numeric(QS)
  if (is.na(QS) || QS < 0) {stop("QS must a positive numeric value")}
  else if (QS < 1 && iside == 1) {stop("For detecting deterioration (side='low') QS must a positive numeric value >= 1")}
  else if (QS > 1 && iside == 2) {stop("For detecting improvement (side='up') QS must a positive numeric value <= 1")}
  itype <- switch(type, cond = 1, cycl = 2)
  if (is.null(itype)) {
    warning("No valid input, using type=cond (conditional steady-state) as default")
    itype <- 1
  }
  m <- as.integer(m)
  if (is.na(m) || m < 0) {stop("m must be a positive integer")}
  .eocusum_adoc_sim(r, k, h, df, coeff, coeff2, QS, iside, itype, m)
}

#' @name eocusum_arloc_h_sim
#' @title Compute alarm threshold of Out of Control EO-CUSUM control charts using simulation
#' @description Compute alarm threshold (Out of Control ARL) of EO-CUSUM control charts using
#'  simulation.

#' @param L0 double. Prespecified in-control Average Run Length.
#' @param df DataFrame. First column are Parsonnet Score values within a range of zero to 100 representing
#' the preoperative patient risk. The second column are binary (0/1) outcome values of each operation.
#' @param k double. Reference value of the CUSUM control chart.
#' @param m integer. Number of simulation runs.
#' @param QS double. Defines the performance of a surgeon with the odds ratio ratio of death Q.
#' @param side character. Default is "low" to calculate ARL for the upper arm of the V-mask. If side = "up",
#'  calculate the lower arm of the V-mask.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model.
#' @param coeff2 NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta}
#'  from the binary logistic regression model of a resampled dataset.
#' @param nc integer. Number of cores.
#' @param verbose boolean. If TRUE verbose output is included, if FALSE a quiet calculation of h is done.
#'
#' @return Returns a single value which is the control limit h for a given ARL.
#'
#' @details The function \code{eocusum_arloc_h_sim} determines the control limit for given in-control ARL (L0) by applying a
#' multi-stage search procedure which includes secant rule and the parallel version of \code{\link{eocusum_arloc_sim}}
#' using \code{\link{mclapply}}.
#'
#' @template eocusum_arloc_h_sim
#'
#' @author Philipp Wittenberg
#' @export
eocusum_arloc_h_sim <- function(L0, k, df, coeff, coeff2, m = 100, QS = 1, side = "low", nc = 1, verbose = FALSE) {
  L0 <- as.integer(L0)
  if (is.na(L0) || L0 <= 0) {stop("Given in-control ARL 'L0' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  if (is.null(coeff2) || is.na(coeff2) || length(coeff2) != 2) {stop("Model coefficients 'coeff2' must be a numeric vector with two elements")}
  coeff2 <- as.vector(coeff2)
  side <- switch(as.character(side), low = 1, up = 2)
  if (is.null(side)) {
    warning("no valid input, using side=low (deterioration) as default")
    side <- 1
  }
  QS <- as.numeric(QS)
  if (is.na(QS) || QS < 0) {stop("QS must a positive numeric value")}
  else if (QS < 1 && side == 1) {stop("For detecting deterioration (side='low') QS must a positive numeric value >= 1")}
  else if (QS > 1 && side == 2) {stop("For detecting improvement (side='up') QS must a positive numeric value <= 1")}
  h2 <- 1
  L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arloc_sim, h = h2, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
  if ( verbose ) cat(paste("(i)\t", h2, "\t", L2, "\n"))
  LL <- NULL
  while ( L2 < L0 & h2 < 6 ) {
    L1 <- L2
    h2 <- h2 + 1
    L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arloc_sim, h = h2, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
    if ( verbose ) cat(paste("(ii)\t", h2, "\t", L2, "\n"))
    LL <- c(LL, L2)
  }
  if ( L2 < L0 ) {
    x  <- 1:5
    LM <- stats::lm(LL ~ I(x) + I(x^2) )
    beta <- stats::coef(LM)
    p <- beta[2] / beta[3]
    q <- (beta[1] - L0)/beta[3]
    h2 <- -p / 2 + 1*sqrt(p^2 / 4 - q)
    L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arloc_sim, h = h2, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
    if ( verbose ) cat(paste("(iii)\t", h2, "\t", L2, "\n"))
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        h2 <- h2 + 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arloc_sim, h = h2, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
        if ( verbose ) cat(paste("(iv)a\t", h2, "\t", L2, "\n"))
        }
      h1 <- h2 - 1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        h2 <- h2 - 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arloc_sim, h = h2, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
        if ( verbose ) cat(paste("(iv)b\t", h2, "\t", L2, "\n"))
      }
      h1 <- h2 + 1
    }
    } else {
    h1 <- h2 - 1
  }
  h.error <- 1
  a.error <- 1
  scaling <- 10^3
  while ( a.error > 1e-4 & h.error > 1e-6 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arloc_sim, h = h3, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
    if ( verbose ) cat(paste("(v)\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    a.error <- min( c(abs(L2 - L0), abs(L2 - L1) ) )
    h.error <- abs(h2 - h1)
    if ( h.error < 0.5 / scaling ) {
      if ( L3 < L0 ) {
        h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
        L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arloc_sim, h = h3, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
        if ( verbose ) cat(paste("(vi)\t", h3, "\t", L3, "\n"))
      }
      break
    }
  }
  if ( L3 < L0 ) {
    h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
    L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arloc_sim, h = h3, k = k, df = df, QS = QS, side = side, coeff = coeff, coeff2 = coeff2, mc.cores = nc)))
    if ( verbose ) cat(paste("(vii)\t", h3, "\t", L3, "\n"))
  }
  h <- h3
  h
}

#' @name eocusum_arl_h_sim
#' @title Compute alarm threshold of EO-CUSUM control charts using simulation
#' @description Compute alarm threshold of EO-CUSUM control charts using simulation.
#'
#' @param L0 double. Prespecified in-control Average Run Length.
#' @param k double. Reference value of the CUSUM control chart.
#' @param m integer. Number of simulation runs.
#' @param side character. Default is "low" to calculate ARL for the upper arm of the V-mask. If side = "up",
#'  calculate the lower arm of the V-mask.
#' @param df DataFrame. First column are Parsonnet Score values within a range of zero to 100 representing
#' the preoperative patient risk. The second column are binary (0/1) outcome values of each operation.
#' @param coeff NumericVector. Estimated coefficients \eqn{\alpha}{alpha} and \eqn{\beta}{beta} from the binary
#' logistic regression model. For more information see details.
#' @param yemp boolean. Use emirical outcome value.
#' @param nc integer. Number of cores.
#' @param verbose boolean. If TRUE verbose output is included, if FALSE a quiet calculation of h is done.
#'
#' @return Returns a single value which is the control limit h for a given ARL.
#'
#' @template eocusum_arl_h_sim
#'
#' @details The function \code{eocusum_arl_h_sim} determines the control limit for given in-control ARL (L0) by applying a
#' multi-stage search procedure which includes secant rule and the parallel version of \code{\link{eocusum_arl_sim}}
#' using \code{\link{mclapply}}.
#'
#' @author Philipp Wittenberg
#' @export
eocusum_arl_h_sim <- function(L0, k, df, coeff, m = 100, yemp = TRUE, side = "low", nc = 1, verbose = FALSE) {
  L0 <- as.integer(L0)
  if (is.na(L0) || L0 <= 0) {stop("Given in-control ARL 'L0' must be a positive integer")}
  k <- as.numeric(k)
  if (is.na(k) || k  < 0) {stop("Reference value 'k' must be a positive numeric value")}
  if (class(df) != "data.frame") {stop("Provide a dataframe for argument 'df'")}
  else if (ncol(df) != 2) {stop("Provide a dataframe with two columns for argument 'df'")}
  else if (vapply(df, class, "")[1] != "integer") {stop("First column of dataframe must be of type integer")}
  else if (vapply(df, class, "")[2] != "numeric") {stop("Second column of dataframe must be of type numeric")}
  df <- as.data.frame(df)
  if (is.null(coeff) || is.na(coeff)  || length(coeff)  != 2) {stop("Model coefficients 'coeff' must be a numeric vector with two elements")}
  coeff <- as.vector(coeff)
  side <- switch(as.character(side), low = 1, up = 2)
  if (is.null(side)) {
    warning("no valid input, using side=low (deterioration) as default")
    side <- 1
  }
  h2 <- 1
  L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, k = k, h = h2, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
  LL <- NULL
  while ( L2 < L0 & h2 < 6 ) {
    L1 <- L2
    h2 <- h2 + 1
    L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, k = k, h = h2, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
    if ( verbose ) cat(paste("(ii)\t", h2, "\t", L2, "\n"))
    LL <- c(LL, L2)
  }
  if ( L2 < L0 ) {
    x  <- 1:5
    LM <- stats::lm(LL ~ I(x) + I(x^2) )
    beta <- stats::coef(LM)
    p <- beta[2] / beta[3]
    q <- (beta[1] - L0) / beta[3]
    h2 <- -p / 2 + 1 * sqrt(p^2 / 4 - q)
    L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, k = k, h = h2, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
    if ( verbose ) cat(paste("(iii)\t", h2, "\t", L2, "\n"))
    if ( L2 < L0 ) {
      while ( L2 < L0 ) {
        L1 <- L2
        h2 <- h2 + 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, k = k, h = h2, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
        if ( verbose ) cat(paste("(iv)a\t", h2, "\t", L2, "\n"))
        }
      h1 <- h2 - 1
    } else {
      while ( L2 >= L0 ) {
        L1 <- L2
        h2 <- h2 - 1
        L2 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, k = k, h = h2, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
        if ( verbose ) cat(paste("(iv)b\t", h2, "\t", L2, "\n"))
      }
      h1 <- h2 + 1
    }
    } else {
    h1 <- h2 - 1
  }
  h.error <- 1
  a.error <- 1
  scaling <- 10^3
  while ( a.error > 1e-4 & h.error > 1e-6 ) {
    h3 <- h1 + (L0 - L1) / (L2 - L1) * (h2 - h1)
    L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, k = k, h = h3, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
    if ( verbose ) cat(paste("(v)\t", h3, "\t", L3, "\n"))
    h1 <- h2
    h2 <- h3
    L1 <- L2
    L2 <- L3
    a.error <- min( c(abs(L2 - L0), abs(L2 - L1) ) )
    h.error <- abs(h2 - h1)
    if ( h.error < 0.5 / scaling ) {
      if ( L3 < L0 ) {
        h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
        L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, k = k, h = h3, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
        if ( verbose ) cat(paste("(vi)\t", h3, "\t", L3, "\n"))
      }
      break
    }
  }
  if ( L3 < L0 ) {
    h3 <- ( round( h3 * scaling ) + 1 ) / scaling - 1e-6
    L3 <- mean(do.call(c, parallel::mclapply(1:m, eocusum_arl_sim, k = k, h = h3, df = df, yemp = yemp, side = side, coeff = coeff, mc.cores = nc)))
    if ( verbose ) cat(paste("(vii)\t", h3, "\t", L3, "\n"))
  }
  h <- h3
  h
}
