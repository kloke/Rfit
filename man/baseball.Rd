\name{baseball}
\alias{baseball}
\docType{data}
\title{Baseball Card Data}
\description{
These data come from the back-side of 59 baseball cards that Carrie had.
}
\usage{data(baseball)}
\format{
  A data frame with 59 observations on the following 6 variables.
  \describe{
    \item{\code{height}}{Height in inches}
    \item{\code{weight}}{Weight in pounds}
    \item{\code{bat}}{a factor with levels \code{L} \code{R} \code{S}}
    \item{\code{throw}}{a factor with levels \code{L} \code{R}}
    \item{\code{field}}{a factor with levels \code{0} \code{1}}
    \item{\code{average}}{ERA if the player is a pitcher and his batting average if the player is a fielder}
  }
}
\source{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\examples{
data(baseball)
wilcox.test(height~field,data=baseball)
rfit(weight~height,data=baseball)
}
\keyword{datasets}
