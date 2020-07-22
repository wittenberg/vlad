#' @references Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455--2473.
#' @examples
#' \dontrun{
#' library("dplyr")
#' library("tidyr")
#' library(ggplot2)
#'
#' ## Datasets
#' data("cardiacsurgery", package = "spcadjust")
#' cardiacsurgery <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0))
#' s5000 <- sample_n(cardiacsurgery, size = 5000, replace = TRUE)
#' df1 <- select(cardiacsurgery, s, y)
#' df2 <- select(s5000, s, y)
#'
#' ## estimate coefficients from logit model
#' coeff1 <- round(coef(glm(y ~ s, data = df1, family = "binomial")), 3)
#' coeff2 <- round(coef(glm(y ~ s, data = df2, family = "binomial")), 3)
#'
#' ## set up
#' RNGkind("L'Ecuyer-CMRG")
#' m <- 10^3
#' kopt <- optimal_k(QA = 2, df = S2I, coeff = coeff1, yemp = FALSE)
#' h <- eocusum_arloc_h_sim(L0 = 370, df = df1, k = kopt, m = m, side = "low", coeff = coeff1,
#'                          coeff2 = coeff2, nc = 4)
#'
#' ## Serial simulation
#' RLS <- do.call(c, lapply(1:m, eocusum_arloc_sim, h = h, k = kopt, df = df1, side = "low",
#'                          coeff = coeff1, coeff2 = coeff2))
#' data.frame(cbind(ARL = mean(RLS), ARLSE = sd(RLS)/sqrt(m)))
#'
#' ## Parallel simulation (FORK)
#' RLS <- simplify2array(parallel::mclapply(1:m, eocusum_arloc_sim, h = h, k = kopt, df = df1,
#'                                          side = "low", coeff = coeff1, coeff2 = coeff2,
#'                                          mc.cores = parallel::detectCores()))
#' data.frame(cbind(ARL = mean(RLS), ARLSE = sd(RLS)/sqrt(m)))
#'
#' ## Parallel simulation (PSOCK)
#' no_cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(no_cores)
#' side <- "low"
#' h_vec <- h
#' QS_vec <- 1
#' k <- kopt
#' parallel::clusterExport(cl, c("h_vec", "eocusum_arloc_sim", "df1", "coeff1", "coeff2",
#'                               "QS_vec", "side", "k"))
#' time <- system.time( {
#'   RLS <- array(NA, dim = c( length(QS_vec), length(h_vec), m))
#'   for (h in h_vec) {
#'     for (QS in QS_vec) {
#'       cat(h, " ", QS, "\n")
#'       RLS[which(QS_vec==QS), which(h==h_vec), ] <- parallel::parSapply(cl, 1:m, eocusum_arloc_sim,
#'                                                                        side = side, QS = QS, h = h,
#'                                                                        k = k, df = df1,
#'                                                                        coeff = coeff1,
#'                                                                        coeff2 = coeff2,
#'                                                                        USE.NAMES = FALSE)
#'     }
#'   }
#' } )
#' ARL <- apply(RLS, c(1, 2), mean)
#' ARLSE <- sqrt(apply(RLS, c(1, 2), var)/m)
#' print(list(ARL, ARLSE, time))
#' parallel::stopCluster(cl)
#' }
