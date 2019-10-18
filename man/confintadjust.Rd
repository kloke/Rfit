\name{confintadjust}
\alias{confintadjust}
\alias{confintadjust.methods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence interval adjustment methods
}
\description{
Returns the critical value to be used in calculating adjusted confidence intervals.
Currently provides methods for Boneferroni and Tukey for confidence interval adjustment methods as well as no adjustment.
}
\usage{
confintadjust(n, k, alpha = 0.05, method = confintadjust.methods, ...)
}
\arguments{
  \item{n}{
sample size
}
  \item{k}{
number of comparisons
}
  \item{alpha}{
overall (experimentwise) type I error rate
}
  \item{method}{
one of confintadjust.methods
}
  \item{\dots}{
Additonal arguments.  Currently not used.
}
}
\details{
Returns critial value based on one of the adjustment methods.
}
\value{
\item{cv}{critical value}
\item{method}{the method used}
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\author{
Joseph McKean, John Kloke
}

\seealso{
\code{\link{oneway.rfit}}
}
