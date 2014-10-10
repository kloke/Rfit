\name{telephone}
\alias{telephone}
\docType{data}
\title{ Telephone Data}
\description{
The number of telephone calls (in tens of millions) made in Belgium from 1950-1973.

}
\usage{data(telephone)}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{\code{year}}{years since 1950 AD}
    \item{\code{calls}}{number of telephone calls in tens of millions}
  }
}
%\details{
  %~~ If necessary, more details than the __description__ above ~~
%}
\source{
Rousseeuw, P.J. and Leroy, A.M. (1987), \emph{Robust Regression and Outlier Detection}, New York: Wiley.
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\examples{
data(telephone)
plot(telephone)
abline(rfit(calls~year,data=telephone))
}
\keyword{datasets}
