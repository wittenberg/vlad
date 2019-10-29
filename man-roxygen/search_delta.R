#' @references Knoth S, Wittenberg P and Gan FF (2019).
#' Risk-adjusted CUSUM charts under model error. To appear in
#' \emph{Statistics in Medicine}, \doi{10.1002/sim.8104}.
#' @examples
#' \dontrun{
#' ## load data
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## preprocess data to 30 day mortality and subset data to
#' ## phase I (In-control) and phase II (monitoring)
#' SALL <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II")))
#'
#' ## subset phase I (In-control)
#' SI <- filter(SALL, phase == "I") %>% select(s, y)
#'
#' ## search delta
#' dML <- search_delta(SI$s, SI$y, type = "ML")
#' dQQ <- search_delta(SI$s, SI$y, type = "Pearson")
#'
#' ## show Log-likelihood (ell()) and Pearson measure (QQ()) for each delta
#' delta <- c(-2, -1, 0, dML, dQQ, 0.5, 1, 2)
#' r <- sapply(delta, function(i) rbind(i, ell(SI$s, SI$y, i), QQ(SI$s, SI$y, i)))
#' rownames(r) <- c("d", "l", "S")
#' t(r)
#' data.frame(t(r)) %>% filter(l == max(l) | S == min(S))
#' }
