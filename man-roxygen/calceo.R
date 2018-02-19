#' @references Lovegrove J, Valencia O, Treasure T, Sherlaw-Johnson C and Gallivan S (1997).
#'  “Monitoring the results of cardiac surgery by variable life-adjusted display.”
#'  \emph{The Lancet}, \strong{350}(9085), pp. 1128–1130.
#'  doi: \href{https://doi.org/10.1016/S0140-6736(97)06507-0}{10.1016/S0140-6736(97)06507-0}.
#'
#' Poloniecki J, Valencia O and Littlejohns P (1998).
#'  “Cumulative risk adjusted mortality chart for detecting changes in death rate: observational
#'  study of heart surgery.” \emph{BMJ}, \strong{316}(7146), pp. 1697–1700.
#'  doi: \href{https://doi.org/10.1136/bmj.316.7146.1697}{10.1136/bmj.316.7146.1697}.
#'
#' Steiner S (2014). “Risk-Adjusted Monitoring of Outcomes in Health Care.”
#'  In Lawless JF (ed.), \emph{Statistics in Action}, pp. 225-242. Informa UK Limited.
#'  doi: \href{https://doi.org/10.1201/b16597-15}{10.1201/b16597-15}.
#'
#' @examples
#' library("vlad")
#' # see Steiner (2014) p. 234
#' coeff <- c("(Intercept)"=-3.68, "Parsonnet"=0.077)
#' # penalty reward for death (E-O scores multiplied with -1 to get O-E scores)
#' calceo(df=data.frame(as.integer(0), 1), coeff=coeff)*-1
#' calceo(df=data.frame(as.integer(50), 1), coeff=coeff)*-1
#' # penalty reward for survival
#' calceo(df=data.frame(as.integer(0), 0), coeff=coeff)*-1
#' calceo(df=data.frame(as.integer(50), 0), coeff=coeff)*-1
#'
#' # Plot a VLAD/CRAM chart
#' library("spcadjust")
#' data("cardiacsurgery")
#' cardiacsurgery <- dplyr::mutate(cardiacsurgery, phase=factor(ifelse(date < 2*365, "I", "II")))
#' S2 <- subset(cardiacsurgery, c(surgeon==2), c("phase", "Parsonnet", "status"))
#' S2I <- subset(S2, c(phase=="I"))
#' S2II <- subset(S2, c(phase=="II"))
#' coeff <- coef(glm(status ~ Parsonnet, data=S2I, family="binomial"))
#' EO <- sapply(1:nrow(S2), function(i) calceo(df=S2[i, c("Parsonnet", "status")], coeff=coeff))
#'
#' df1 <- data.frame(cbind(subset(S2, select=c("phase")), "n"=1:nrow(S2), "cEO"=cumsum(EO)))
#' df2 <- tidyr::gather(df1, "variable", value, c(-n, -phase))
#'
#' ggplot2::qplot(data=df2, n, value, colour=phase, geom=c("line", "point"),
#'                xlab="Patient number", ylab="CUSUM E-O") +
#'   ggplot2::geom_hline(yintercept=0, linetype="dashed") +
#'   ggplot2::theme_classic()
