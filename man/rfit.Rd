\name{rfit}
\alias{rfit}
\alias{rfit.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rank-based Estimates of Regression Coefficients }
\description{
Minimizes Jaeckel's dispersion function to obtain a rank-based solution for linear models.
}
\usage{
rfit(formula, data = list(), ...)
%rfit.formula(formula, data = list(), ...)
\method{rfit}{default}(formula, data, subset, yhat0 = NULL, 
scores = Rfit::wscores, symmetric = FALSE, TAU = "F0", 
betahat0 = NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class formula }
  \item{data}{ an optional data frame}
  \item{subset}{ an optional argument specifying the subset of observations to be used}
  \item{yhat0}{ an n by 1 vector of initial fitted values, default is NULL}
  \item{scores}{ an object of class 'scores'}
  \item{symmetric}{ logical.  If 'FALSE' uses median of residuals as estimate of intercept} 
  \item{TAU}{version of estimation routine for scale parameter.  F0 for Fortran, R for (slower) R, N for none}
  \item{betahat0}{ a p by 1 vector of initial parameter estimates, default is NULL}
  \item{\dots}{ additional arguments to be passed to fitting routines }

}
\details{
Rank-based estimation involves replacing the L2 norm of least squares estimation with a pseudo-norm which is a function of the residuals and the scored ranks of the residuals.
That is, in rank-based estimation, the usual notion of Euclidean distance is replaced with another measure of distance which is referred to as Jaeckel's (1972) dispersion function.
Jaeckel's dispersion function depends on a score function and a library of commonly used score functions is included; eg., linear (Wilcoxon) and normal (Gaussian) scores.
If an inital fit is not supplied (i.e. yhat0 = NULL and betahat0 = NULL) then inital fit is based on a LS fit.

Esimation of scale parameter tau is provided which may be used for inference.  
}
\value{
\item{coefficients}{estimated regression coefficents with intercept}
\item{residuals}{the residuals, i.e. y-yhat}
\item{fitted.values}{ yhat = x betahat}
%\item{scores}{ score function used in estimation}
%\item{x}{original design matrix}
%\item{y}{original response vector}
\item{xc}{centered design matrix}
\item{tauhat}{estimated value of the scale parameter tau}
\item{taushat}{estimated value of the scale parameter tau_s}
%\item{symmetric}{
\item{betahat}{estimated regression coefficents}
\item{call}{Call to the function}
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Jaeckel, L. A. (1972). Estimating regression coefficients by minimizing the dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449 - 1458.

Jureckova, J. (1971). Nonparametric estimate of regression coefficients. \emph{Annals of Mathematical Statistics}, 42, 1328 - 1338.
}
\author{ John Kloke, Joesph McKean }
\seealso{
\code{\link{summary.rfit}}
\code{\link{drop.test}}
\code{\link{rstudent.rfit}}
}

\examples{
data(baseball)
data(wscores)
fit<-rfit(weight~height,data=baseball)
summary(fit)

### set the starting value
x1 <- runif(47); x2 <- runif(47); y <- 1 + 0.5*x1 + rnorm(47)
# based on a fit to a sub-model
rfit(y~x1+x2,yhat0=fitted.values(rfit(y~x1)))

### set value of delta used in estimation of tau ###
w <- factor(rep(1:3,each=3))
y <- rt(9,9)
rfit(y~w)$tauhat
rfit(y~w,delta=0.95)$tauhat  # recommended when n/p < 5
}

\keyword{nonparametric}
\keyword{robust}
\keyword{regression}

