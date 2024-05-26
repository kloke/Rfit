\name{summary.rfit}
\alias{summary.rfit}
\title{ Summarize Rank-Based Linear Model Fits }
\description{
Provides a summary similar to the traditional least squares fit.
}
\usage{
\method{summary}{rfit}(object,overall.test,...)
}
\arguments{
\item{object}{an object of class 'rfit', usually, a result of a call to 'rfit'}
\item{overall.test}{either 'wald' or 'drop'}
\item{\dots}{ additional arguments }
}
\details{
Provides summary statistics based on a rank-based fit.
A table of estimates, standard errors, t-ratios, and p-values are provided.
An overall test of the explantory variables is provided; the default is to use a Wald test.
A drop in dispersion test is also availble in which case a robust R^2 is provided as well.
}
%\value{
  %~Describe the value returned
  %If it is a LIST, use
  %\item{comp1 }{Description of 'comp1'}
  %\item{comp2 }{Description of 'comp2'}
  %...
%}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\author{ John Kloke }


% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(baseball)
fit<-rfit(weight~height,data=baseball)
summary(fit)
summary(fit,overall.test='drop')
}
