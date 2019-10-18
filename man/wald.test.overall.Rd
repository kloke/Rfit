\name{wald.test.overall}
\alias{wald.test.overall}
\title{Overall Wald test}
\description{Conducts a Wald test of all regression parameters are zero}
\usage{wald.test.overall(fit)}
\arguments{
\item{fit}{result from a rfit}
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\author{John Kloke}
\examples{
x <- rnorm(47)
y <- rnorm(47)
wald.test.overall(rfit(y~x))
}
\keyword{htest}
\keyword{robust}
\keyword{regression}
