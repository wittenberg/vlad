#' @references Wittenberg P, Gan FF, Knoth S (2018).
#' A simple signaling rule for variable life-adjusted display derived from
#' an equivalent risk-adjusted CUSUM chart.
#' \emph{Statistics in Medicine}, \strong{37}(16), pp 2455-2473.
#'
#' @examples
#' \dontrun{
#' library("dplyr")
#' library("tidyr")
#' library(ggplot2)
#' data("cardiacsurgery", package = "spcadjust")
#'
#' ## preprocess data to 30 day mortality and subset phase I (In-control) of surgeons 2
#' SALL <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'          phase = factor(ifelse(date < 2*365, "I", "II")))
#'
#' ## subset phase I (In-control)
#' SI <- filter(SALL, phase == "I") %>% select(s, y)
#'
#' ## estimate coefficients from logit model
#' coeff1 <- coef(glm(y ~ s, data = SI, family = "binomial"))
#'
#' ## determine k for detecting deterioration
#' kopt <- optimal_k(QA = 2, df = SI, coeff = coeff, yemp = FALSE)
#'
#' ## subset phase II of surgeons 2
#' S2II <- filter(SALL, phase == "II", surgeon == 2) %>% select(s, y)
#' n <- nrow(S2II)
#'
#' ## CUSUM statistic without reset
#' z <- sapply(1:n, function(i) calceo(df = S2II[i, c("s", "y")], coeff = coeff1))
#' cv <- eocusum_scores(z = z, k = kopt)
#' s1 <- cv$s1; s1l <- cv$s1l
#' dm1 <- data.frame(cbind("n" = 1:length(s1), "Cup" = s1, "Clow" = s1l, "h1" = 2, "h2" = -2))
#'
#' ## CUSUM statistic reset after signal
#' cv <- eocusum_scores(z = z, k = kopt, reset = TRUE, h1 = 2, h2 = 2)
#' s1 <- cv$s1; s1l <- cv$s1l
#' dm2 <- data.frame(cbind("n" = 1:length(s1), "Cup" = s1, "Clow" = s1l, "h1" = 2, "h2" = -2))
#'
#' dm3 <- bind_rows(dm1, dm2, .id = "type")
#' dm3$type <- recode_factor(dm3$type, `1`="No resetting", `2`="Resetting")
#' dm3 %>%
#'   gather("CUSUM", value, c(-n, - type)) %>%
#'   ggplot(aes(x = n, y = value, colour = CUSUM, group = CUSUM)) +
#'   geom_hline(yintercept = 0, colour = "darkgreen", linetype = "dashed") +
#'   geom_line(size = 0.5) +
#'   facet_wrap( ~ type, ncol = 1, scales = "free") +
#'   labs(x = "Patient number n", y = "CUSUM values") + theme_classic() +
#'   scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
#'   scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
#'   guides(colour = "none") +
#'   scale_color_manual(values = c("blue", "orange", "red", "red"))
#' }
