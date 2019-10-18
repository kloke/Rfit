\name{rstudent.rfit}
\alias{rstudent.rfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Studentized Residuals for Rank-Based Regression }
\description{ Returns the Studentized residuals based on rank-based estimation.  }
\usage{
\method{rstudent}{rfit}(model,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ an object of class rfit }
  \item{\dots}{ additional arguments. currently not used.}
}
%\details{
  %~~ If necessary, more details than the description above ~~
%}
%\value{
  %~Describe the value returned
  %If it is a LIST, use
  %\item{comp1 }{Description of 'comp1'}
  %\item{comp2 }{Description of 'comp2'}
  %...
%}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }

\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}

\author{ John Kloke, Joseph McKean}
\seealso{ \code{\link{rfit}}}
\examples{
x<-runif(47)
y<-rcauchy(47)
qqnorm(rstudent(fit<-rfit(y~x)))
plot(x,rstudent(fit)) ; abline(h=c(-2,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{robust}
