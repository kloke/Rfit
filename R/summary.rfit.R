#' Summarize Rank-Based Linear Model Fits
#' 
#' Provides a summary similar to the traditional least squares fit.
#' 
#' 
#' @param object an object of class 'rfit', usually, a result of a call to
#' 'rfit'
#' @param \dots additional arguments
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @examples
#' 
#' data(baseball)
#' fit<-rfit(weight~height,data=baseball)
#' summary(fit)
#' 
#' @export summary.rfit
summary.rfit <- function (object,overall.test='all',...) {

  tauhat <- object$tauhat
  n<-length(object$y)
  pp1 <- object$qrx1$rank
  est <- object$coef
  ses <- sqrt(diag(vcov(object)))
  tstat <- est/ses
  pval <- 2 * pt(-abs(tstat), n - pp1)
  coef <- cbind(est, ses, tstat, pval)
  colnames(coef) <- c("Estimate", "Std. Error", "t.value","p.value")

  ans <- list(coefficients = coef)

  if( overall.test == 'all' || overall.test == 'wald' ) {
    wt <- wald.test.overall(object)
    ans <- utils::modifyList(ans,list(waldstat = wt$F, waldpval = wt$p.value))
  } 
  if( overall.test == 'all' ||  overall.test == 'drop') {
    dt <- drop.test(object)
    R2 <- (dt$df1/dt$df2 * dt$F)/(1 + dt$df1/dt$df2 * dt$F)
    ans <- utils::modifyList(ans,list(dropstat = dt$F, droppval = dt$p.value, R2 = R2))
  }
  ans$overall.test <- overall.test
  ans$call <- object$call
  class(ans) <- "summary.rfit"
  ans

}
