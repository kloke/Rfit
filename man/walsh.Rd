\name{walsh}
\alias{walsh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Walsh Averages }
\description{
Given a list of n numbers, the Walsh averages are the \eqn{latex}{ n(n+1)/2 } pairwise averages.
}
\usage{
walsh(x)
}
\arguments{
  \item{x}{ A numeric vector  }
}
\value{
The Walsh averages.
}
\references{ 
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Hollander, M. and Wolfe, D.A. (1999), \emph{Nonparametric Statistical Methods}, New York: Wiley.
}
%\note{ ~~further notes~~ 

\author{ John Kloke, Joseph McKean }

 %~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{signedrank}}}
\examples{

median(walsh(rnorm(100)))  # Hodges-Lehmann estimate of location

## The function is currently defined as
function (x) 
{
    n <- length(x)
    w <- vector(n * (n + 1)/2, mode = "numeric")
    ind <- 0
    for (i in 1:n) {
        for (j in i:n) {
            ind <- ind + 1
            w[ind] <- 0.5 * (x[i] + x[j])
        }
    }
    return(w)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{n }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
