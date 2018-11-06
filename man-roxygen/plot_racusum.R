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
#' ## show data
#' head(SALL)
#'
#' ## subset phase I (In-control)
#' SALLI <- filter(SALL, phase == "I") %>% select(s, y)
#'
#' ## estimate coefficients from logit model
#' mod1 <- glm(y ~ s, data = SALLI, family = "binomial")
#' summary(mod1)
#'
#' ## extract coefficients from logit model
#' coeff <- round(coef(mod1), 3)
#'
#' ## subset phase II of surgeons 2
#' Si <- filter(SALL, surgeon == 2, phase == "II") %>% select(s, y)
#'
#' ## compute control limit for detecting deterioration RA = 2 / improvement RA = 1/2
#' ## setup
#' L0 <- 200
#' m <- 10^4
#' nc <- 4
#' UCL <- racusum_crit_sim(L0 = L0, df = SALLI, coeff = coeff, m = m, nc = nc, RA = 2)
#' LCL <- racusum_crit_sim(L0 = L0, df = SALLI, coeff = coeff, m = m, nc = nc, RA = 1/2)
#' print(c(UCL, LCL))
#'
#' ## plot RA-CUSUM (LLR) control chart with reset and signals
#' p1 <- plot_racusum(data = Si, coeff = coeff, reset = TRUE,  h1 = UCL, h2 = LCL, signal = TRUE)
#' p1 + guides(colour = "none")
#'
#' ## show alarms
#' filter(p1$data, signal %in% c("AlarmU", "AlarmL"))
#' ## plot RA-CUSUM (LLR) control chart without reset and signals
#' p2 <- plot_racusum(data = Si, coeff = coeff, reset = TRUE,  h1 = UCL, h2 = LCL, signal = FALSE)
#' p2
#'
#' ## show first signals
#' filter(p2$data, signal %in% c("AlarmU", "AlarmL")) %>% head()
#' }
