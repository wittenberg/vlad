#' @references Steiner SH, Cook RJ, Farewell VT and Treasure T (2000).
#' Monitoring surgical performance using risk-adjusted cumulative sum charts.
#'  \emph{Biostatistics}, \strong{1}(4), pp. 441--452.
#'
#' Steiner S (2014). Risk-Adjusted Monitoring of Outcomes in Health Care.
#' In Lawless JF (ed.), \emph{Statistics in Action}, pp. 225--242. Informa UK Limited.
#'
#' Parsonnet V, Dean D, Bernstein AD (1989). A method of uniform stratification of risk
#' for evaluating the results of surgery in acquired adult heart disease.
#' \emph{Circulation}, \strong{79}(6):I3--12.
#'
#' Rigdon SE and Fricker RD (2015). Health Surveillance.
#' In Chen DG and Wilson J (eds) \emph{Innovative Statistical Methods for Public Health Data},
#' pp. 203--249. Springer, Cham.
#'
#' @examples
#' \dontrun{
#' library(vlad)
#' ## see Steiner et al. 2000 p. 445 or Steiner (2014) p. 234
#' coeff <- c("(Intercept)" = -3.68, "Parsonnet" = 0.077)
#' ## low risk patient (Parsonnet score=0) has a risk of death 2.5%
#' gettherisk(0L, coeff = coeff)
#' ## high risk patient (Parsonnet score=71) has a risk of death 86%
#' gettherisk(71L, coeff = coeff)
#' ## high risk patient (Parsonnet score=50) has a risk of death 54%
#' gettherisk(50L, coeff = coeff)
#'
#' ## see Rigdon and Fricker (2015) p. 221 and p. 225
#' coeff <- c("(Intercept)" = -3.67, "Parsonnet" = 0.077)
#' ## patients probability of death 0.09912 for Parsonnet score 19
#' round(gettherisk(19L, coeff), 5)
#' ## patients probability of death 0.02484 for Parsonnet score 0
#' round(gettherisk(0L, coeff), 5)
#'
#' ## preprocess data to 30 day mortality and subset phase I (In-control)
#' library("dplyr")
#' data("cardiacsurgery", package = "spcadjust")
#' SI <- cardiacsurgery %>% rename(s = Parsonnet) %>%
#'   mutate(y = ifelse(status == 1 & time <= 30, 1, 0),
#'         phase = factor(ifelse(date < 2*365, "I", "II"))) %>%
#'   filter(phase == "I") %>% select(s, y)
#'
#' ## Get mortality and probability of death of a phase I dataset
#' GLM1 <- glm(y ~ s, data = SI, family = "binomial")
#' coeff1 <- coef(GLM1)
#' mprob <- as.numeric(table(SI$s) / length(SI$s))
#'
#' ## Use estimated model coefficients and parsonnet scores in function gettherisk()
#' ## or predicted values from a GLM
#' usi <- sort(unique(SI$s))
#' mort <- sapply(usi, gettherisk, coeff = coeff1)
#' mort1 <- predict(GLM1, newdata = data.frame(s = usi), type = "response")
#' all.equal(as.numeric(mort), as.numeric(mort1))
#' df1 <- data.frame(s = usi, mprob, mort)
#'
#' ## Plot mortality and estimated probability to die of phase I data
#' library(ggplot2)
#' qplot(data = df1, s, mprob, geom = c("line", "point")) + theme_classic()
#' xx <- tapply(SI$y, SI$s, sum)
#' nn <- tapply(SI$y, SI$s, length)
#' ll <- binom::binom.confint(xx, nn, conf.level = 0.99, methods = "exact")$lower
#' uu <- binom::binom.confint(xx, nn, conf.level = 0.99, methods = "exact")$upper
#' ybar <- tapply(SI$y, SI$s, mean)
#' ggplot(data = df1, aes(s, mort)) +
#'  geom_point(data = data.frame(s = usi, ybar), aes(s, ybar), inherit.aes = FALSE) +
#'  geom_errorbar(aes(ymax = uu, ymin = ll), width = 0.9, position = "dodge", alpha = 0.3) +
#'  geom_line(colour = "red") + labs(x = "Parsonnet score", y = "Estimated Probability to die") +
#'  theme_classic()
#' }
