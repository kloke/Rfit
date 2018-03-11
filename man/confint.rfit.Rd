\name{confint.rfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence Interval Method for Objects of Class Rfit}
\description{
Implements a method for calculating confidence intervals for robust
linear models of class Rfit. Based on Kloke and McKean (2012:58) and the general design of \code{confint.lm}.
}
\usage{
\method{confint}{rfit}(object, level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{rfit} }
  \item{level}{The confidence level for the interval of interest. Defaults to 0.95}
}
\author{ Gustavo A. Ballen \email{gaballench@gmail.com}}
\seealso{ 
  \code{\link{rfit}}
}

