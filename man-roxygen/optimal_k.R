#' @examples
#' library("vlad"); library("spcadjust")
#' data("cardiacsurgery")
#' cardiacsurgery <- dplyr::mutate(cardiacsurgery, phase=factor(ifelse(date < 2*365, "I", "II")))
#' S2I <- subset(cardiacsurgery, c(surgeon==2 & phase=="I"), c("Parsonnet", "status"))
#' coeff <- coef(glm(status ~ Parsonnet, data=S2I, family="binomial"))
#' kopt <- optimal_k(QA=2, df=S2I, coeff=coeff, yemp=FALSE)
#' kopt ## (Deterioration)
## k_opt = 0.04059649
#' # manually find optimal k for detecting improvement
#' QA <- 2
#' pbar <- mean(sapply(S2I[, 1], gettherisk, coef=coeff))
#' kopt <- pbar * ( QA - 1 - log(QA) ) / log(QA)
#' all.equal(kopt, optimal_k(QA=2, df=S2I, coeff=coeff, yemp=FALSE) )
#' kopt <- optimal_k(QA=1/2, df=S2I, coeff=coeff, yemp=FALSE)
#' kopt ##(Improvement)
#' ## k_opt = 0.02555328
#' ### k = kopt
#' QA <- 1/2
#' # manually find optimal k for detecting improvement
#' pbar <- mean(sapply(S2I[, 1], gettherisk, coef=coeff))
#' kopt <- pbar * ( 1 - QA + log(QA) ) / log(QA)
#' all.equal(kopt, optimal_k(QA=1/2, df=S2I, coeff=coeff, yemp=FALSE) )
