\name{gettau}
\alias{gettau}
\alias{gettauF0}
\title{ Estimate of the scale parameter tau}
\description{
An estimate of the scale parameter tau is needed for the standard errors of the coefficents in rank-based regression.
}
\usage{
gettau(ehat, p, scores = Rfit::wscores, delta = 0.8, hparm = 2, ...)
}
\arguments{
  \item{ehat}{ full model residuals }
  \item{p}{ number of regression coefficents }
  \item{scores}{ object of class scores, defaults to Wilcoxon scores}
  \item{delta}{ confidence level }
  \item{hparm}{ Joe's hparm }
  \item{\dots}{ additional arguments. currently unused }
}
\details{
This is the confidence interval type estimate of the scale parameter tau developed my Koul, Sievers, and McKean (1987).
This estimate is also discussed in Section 3.7.1 of Hettmansperger and McKean (1998).
One of these function is called in rfit.  The default is to use the faster FORTRAN version.
The R version can be more precise in small samples, but also can be much slower especially when sample sizes are large.
}
\value{
Length one numeric object.
}
\references{ 
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Koul, H.L., Sievers, G.L., and McKean, J.W. (1987) An esimator of the scale parameter for the rank analysis of linear models under general score functions, \emph{Scandinavian Journal of Statistics}, 14, 131-141.
}
\author{ Joseph McKean, John Kloke}
\seealso{ \code{\link{rfit}} }
