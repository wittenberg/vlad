#' @examples
#' \dontrun{
#' library(vlad)
#' m <- 1e3
#' RLS <- sapply(1:m, racusum_discretebeta_arl_sim, shape1=1, shape2=3, coeff=c(-3.6798, 0.0768),
#' h=4.5, RA=2,  rs=71+1, RQ=1)
#' data.frame(cbind(ARL=mean(RLS), ARLSE=sd(RLS)/sqrt(m)))
#' }
#'
