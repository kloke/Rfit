\name{quail}
\alias{quail}
\docType{data}
\title{ Quail Data }
\description{
Thirty-nine quail were randomized to one of for treatments for lowering cholesterol.
}
\usage{data(quail)}
\format{
  A data frame with 39 observations on the following 2 variables.
  \describe{
    \item{\code{treat}}{a factor with levels \code{1} \code{2} \code{3} \code{4}}
    \item{\code{ldl}}{a numeric vector}
  }
}
\source{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}

\examples{
data(quail)
boxplot(ldl~treat,data=quail)
}
\keyword{datasets}
