\name{drop.test}
\alias{drop.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Drop (Reduction) in Dispersion Test }
\description{
Given two full model fits, this function performs a reduction in dispersion test.
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
Rank-based inference procedure analogous to the traditional (LS) reduced model test.

The full and reduced model dispersions are calculated.
The reduction in dispersion test, or drop test for short, has an asymptotic chi-sq
distribution.
Simulation studies suggest using F critical values.
The p-value returned is based on a F-distribution with df1 and df2 degrees of freedom
where df1 is the difference in the number of parameters in the fits of fitF and fitR
and df2 is the residual degrees of freedom in the fit fitF.

Both fits are based on a minimization routine.
It is possible that resulting solutions are such that the fitF$disp > fitRdisp.
We recommend starting the full model at the reduced model fit as a way to avoid this situation.
See examples.
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

## try starting the full model at the reduced model fit ##
fitF<-rfit(y~x1+x2,yhat0=fitR$fitted)
drop.test(fitF,fitR)
}
