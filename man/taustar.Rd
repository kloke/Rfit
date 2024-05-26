\name{taustar}
\alias{taustar}
\title{ Estimate of the Scale Parameter taustar }
\description{
An estimate of the scale parameter taustar =  1/(2*f(0)) is needed for the standard error of the intercept in rank-based regression.
}
\usage{
taustar(e, p, conf = 0.95)
}
\arguments{
  \item{e}{n x 1 vector of full model residuals }
  \item{p}{is the number of regression coefficients (without the intercept) }
  \item{conf}{ confidence level of CI used }
}
\details{
Confidence interval estimate of taustar.
See, for example, Hettmansperger and McKean (1998) p.7-8 and p.25-26.
}
\value{
Length-one numeric object containing the estimated scale parameter taustar.
}

\references{ 
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}

\author{ Joseph McKean, John Kloke}
%\note{ ~~further notes~~ 

 %~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{rfit}} }
\examples{
##  This is an internal function.  See rfit for user-level examples.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
