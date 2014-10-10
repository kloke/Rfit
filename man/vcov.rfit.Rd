\name{vcov.rfit}
\alias{vcov.rfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Variance-Covariance Matrix for Rank-Based Regression }
\description{
Returns the variance-covariance matrix of the regression estimates from an object of type rfit.
}
\usage{
\method{vcov}{rfit}(object, intercept = NULL,...)
}
\arguments{
  \item{object}{ an object of type rfit }
  \item{intercept}{ logical. If TRUE include the variance-covariance estimates corresponding to the intercept }
  \item{\dots}{ additional arguments }

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
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}

\author{ John Kloke \email{kloke@biostat.wisc.edu}}

\seealso{\code{\link{rfit}}}

