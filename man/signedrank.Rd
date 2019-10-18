\name{signedrank}
\alias{signedrank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Signed-Rank Estimate of Location (Intercept) }
\description{
Returns the signed-rank estimate of intercept with is equivalent to the Hodges-Lehmann estimate of the residuals.
}
\usage{
signedrank(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ numeric vector } 
}
\value{
  Returns the median of the Walsh averages.
}
\references{ 
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Hollander, M. and Wolfe, D.A. (1999), \emph{Nonparametric Statistical Methods}, New York: Wiley.

}
\author{ John Kloke, Joseph McKean }

\seealso{ \code{\link{walsh}}}
\examples{

## The function is currently defined as
function (x) 
median(walsh(x))
}
