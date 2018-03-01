#' @examples
#' \dontrun{
#' library("vlad")
#' library("spcadjust")
#' data("cardiacsurgery")
#' cardiacsurgery <- dplyr::mutate(cardiacsurgery, phase=factor(ifelse(date < 2*365, "I", "II")))
#' S2 <- subset(cardiacsurgery, c(surgeon==2), c("phase", "Parsonnet", "status"))
#' # subset phase I (In-control) of surgeons 2
#' S2I <- subset(S2, c(phase=="I"), c("Parsonnet", "status"))
#' # estimate coefficients from logit model
#' coeff1 <- round(coef(glm(status~Parsonnet, data=S2I, family="binomial")), 3)
#'
#' racusum_arl_h_sim(L0=740, df=S2I, coeff=coeff1, m=10^2, nc=4)
#'}
