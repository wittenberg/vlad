#' @examples
#' \dontrun{
#' library(vlad)
#' m <- 1e3
#' RLS <- sapply(1:m, racusum_beta_arl_sim, h=4.5, shape1=1, shape2=3, coeff=c(-3.6798, 0.0768),
#' RA = 2, RQ = 1)
#' data.frame(cbind(ARL=mean(RLS), ARLSE=sd(RLS)/sqrt(m)))
#' }
#'
