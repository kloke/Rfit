\name{drop.test}
\alias{drop.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Drop (Reduction) in Dispersion Test }
\description{
Given two full model fits, this function performs a reduction in disperion test.
}
\usage{
drop.test(fitF, fitR = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitF}{ An object of class rfit.  The full model fit. }
  \item{fitR}{ An object of class rfit.  The reduced model fit. }
}
\details{
Rank-based inference proceedure analogous to the traditional (LS) reduced model test.
}
\value{
%  ~Describe the value returned
%  If it is a LIST, use
\item{F}{Value of the F test statistic}
\item{p.value}{The observed significance level of the test (using an F quantile)}
\item{RD}{Reduced model dispersion minus Full model dispersion}
\item{tauhat}{Estimate of the scale parameter (using the full model residuals)}
\item{df1}{numerator degrees of freedom}
\item{df2}{denominator degrees of freedom}
  %\item{comp1 }{Description of 'comp1'}
  %\item{comp2 }{Description of 'comp2'}
%  ...
}
\author{ John Kloke \email{kloke@biostat.wisc.edu}}

\references{ 
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}

\seealso{ \code{\link{rfit}}}
\examples{
y<-rnorm(47)
x1<-rnorm(47)
x2<-rnorm(47)
fitF<-rfit(y~x1+x2)
fitR<-rfit(y~x1)
drop.test(fitF,fitR)
}
