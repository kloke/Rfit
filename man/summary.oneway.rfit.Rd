\name{summary.oneway.rfit}
\alias{summary.oneway.rfit}
\title{
Provides a summary for the oneway anova based on an R fit.
}
\description{
Provides a summary for the oneway anova based on an R fit including a test for main effects as tests for pairwise comparisons.
}
\usage{
\method{summary}{oneway.rfit}(object, alpha=0.05,method=confintadjust.methods,...)
}
\arguments{
\item{object}{an object of class 'oneway.rfit', usually, a result of a call to 'oneway.rfit'}
\item{alpha}{Experimentwise Error Rate}
\item{method}{method used in confidence interval adjustment}
\item{\dots}{ additional arguments }
}

\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}

\author{ John Kloke, Joseph McKean}

\examples{
data(quail)
oneway.rfit(quail$ldl,quail$treat)
}

