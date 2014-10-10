\name{print.rfit}
\alias{print.rfit}
\alias{print.summary.rfit}
\alias{print.drop.test}
\alias{print.oneway.rfit}
\alias{print.summary.oneway.rfit}
\alias{print.raov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rfit Internal Print Functions}
\description{
These functions print the output in a user-friendly manner using the internal R function \code{print}.
}
\usage{
\method{print}{rfit}(x, ...)
\method{print}{summary.rfit}(x, digits = max(5, .Options$digits - 2), ...)
\method{print}{drop.test}(x, digits = max(5, .Options$digits - 2), ...)
\method{print}{oneway.rfit}(x, digits = max(5, .Options$digits - 2), ...)
\method{print}{summary.oneway.rfit}(x, digits = max(5, .Options$digits - 2), ...)
\method{print}{raov}(x, digits = max(5, .Options$digits - 2), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object to be printed }
  \item{digits}{number of digits to display}
  \item{\dots}{  additional arguments to be passed to \code{print} }
}
\author{ John Kloke \email{kloke@biostat.wisc.edu}}
\seealso{ 
\code{\link{rfit}}, \code{\link{summary.rfit}}, \code{\link{drop.test}}
}

