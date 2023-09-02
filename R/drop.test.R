#' Drop (Reduction) in Dispersion Test
#' 
#' Given two full model fits, this function performs a reduction in disperion
#' test.  Given one fit, returns the test comparing to the null model.
#' 
#' Rank-based inference proceedure analogous to the traditional (LS) reduced
#' model test.
#' 
#' @param fitF An object of class rfit.  The full model fit.
#' @param fitR An object of class rfit.  The reduced model fit.
#' @return % ~Describe the value returned % If it is a LIST, use \item{F}{Value
#' of the F test statistic} \item{p.value}{The observed significance level of
#' the test (using an F quantile)} \item{RD}{Reduced model dispersion minus
#' Full model dispersion} \item{tauhat}{Estimate of the scale parameter (using
#' the full model residuals)} \item{df1}{numerator degrees of freedom}
#' \item{df2}{denominator degrees of freedom} %\item{comp1 }{Description of
#' 'comp1'} %\item{comp2 }{Description of 'comp2'} % ...
#' @author John Kloke 
#' @seealso \code{\link{rfit}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @examples
#' 
#' y<-rnorm(47)
#' x1<-rnorm(47)
#' x2<-rnorm(47)
#' fitF<-rfit(y~x1+x2)
#' fitR<-rfit(y~x1)
#' drop.test(fitF,fitR)
#' 
#' @export drop.test
drop.test <- function (fitF, fitR = NULL){

  EPS<-.Machine$double.eps*100  # abs(RD) smaller than this considered 0

  pp1 <- fitF$qrx1$rank

  if (is.null(fitR)) {
#    rd <- disp(rep(0, ncol(fitF$x)), fitF$x, fitF$y, fitF$scores) - 
#      disp(fitF$betahat, fitF$x, fitF$y, fitF$scores)
    rd <- fitF$D0 - fitF$D1
    df1 <- pp1 - 1
  } else {

    # check if reduced model design is a subspace of column space of full model 
    Q_R <- qr.Q(fitR$qrx1)[,1:fitR$qrx1$rank]
    if( !all(abs( qr.fitted(fitF$qrx1,Q_R) - Q_R ) < .Machine$double.eps ^ 0.5 ) ) stop('Not apparent reduced model is a proper subset of full model.  Check model assumptions.')
#    if( !all(abs( qr.fitted(fitF$qrx1,qr.Q(fitR$qrx1)) - qr.Q(fitR$qrx1) ) < .Machine$double.eps ^ 0.5 ) ) warning('Not apparent reduced model is a proper subset of full model.  Check model assumptions.')
#    rd <- disp(fitR$betahat, fitR$x, fitR$y, fitR$scores) - 
#      disp(fitF$betahat, fitF$x, fitF$y, fitF$scores)

    rd <- fitR$D1 - fitF$D1
    df1 <- length(fitF$betahat) - length(fitR$betahat)
    if(df1 < 1) stop('Not positive integer value of numerator degrees of freedom found.  Check function call or reduced model assumptions.')
  }

  if( abs(rd) < EPS ) rd <- 0

  if( rd < 0 ) stop( "drop.test: negative reduction in dispersion found\n",
	"try starting full model at reduced model\n",
	"see help(drop.test) for more information" )

  df2 <- length(fitF$y) - pp1
  test <- (rd/df1)/(fitF$tauhat/2)
#  if( abs(test) < EPS ) test <- 0
  pval <- 1 - pf(test, df1, df2)
  ans <- list(F = test, p.value = pval, RD = rd, tauhat = fitF$tauhat, 
    df1 = df1, df2 = df2)
  class(ans) <- "drop.test"
  ans

}
