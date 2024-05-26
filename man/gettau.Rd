\name{gettau}
\alias{gettau}
\alias{gettauF0}
\title{Estimate of the scale parameter tau}
\description{
An estimate of the scale parameter tau may be used for the standard errors of the coefficients in rank-based regression.
}
\usage{
gettau(ehat, p, scores = Rfit::wscores, delta = 0.8, hparm = 2, ...)
gettauF0(ehat, p, scores = Rfit::wscores, delta = 0.8, hparm = 2, ...)
}
\arguments{
  \item{ehat}{vector of length n: full model residuals }
  \item{p}{scalar: number of regression coefficients (excluding the intercept); see Details}
  \item{scores}{ object of class scores, defaults to Wilcoxon scores}
  \item{delta}{ confidence level; see Details}
  \item{hparm}{ used in Huber's degrees of freedom correction; see Details}
  \item{\dots}{ additional arguments. currently unused }
}
\details{
For rank-based analyses of linear models, the estimator \eqn{\hat{\tau}} of the scale parameter \eqn{\tau} plays a standardizing role in the standard errors (SE) of the rank-based estimators of the regression coefficients and in the denominator of Wald-type and the drop-in-dispersion test statistics of linear hypotheses.
\code{rfit} currently implements the KSM (Koul, Sievers, and McKean 1987) estimator of tau.   

The functions \code{gettau} and \code{gettauF0} are both available to compute the KSM estimate and may be call from \code{rfit} and used for inference.  The default is to use the faster FORTRAN version \code{gettauF0} via the to option \code{TAU='F0'}.
The R version, \code{gettau}, may be much slower especially when sample sizes are large; this version may be called from \code{rfit} using the option \code{TAU='R'}.

The KSM estimator tauhat is a density type estimator that has the bandwidth given by \eqn{t_\delta/sqrt{n}},
where \eqn{t_\delta} is the \eqn{\delta-th} quantile of the cdf \eqn{H(y)} given in expression (3.7.2) of Hettmansperger and McKean (2011), with the corresponding estimator \eqn{\hat{H}}, given in expression (3.7.7) of Hettmansperger and McKean (2011).

Based on simulation studies, most situations where (n/p >= 6), the default delta = 0.80 provides a valid rank-based
analysis (McKean and Sheather, 1991).
For situations with n/p < 6, caution is needed as the KSM estimate is sensitive to choice of bandwidth.  McKean and Sheather (1991) recommend using a value of 0.95 for delta in such situations.

To correct for heavy-tailed random errors, Huber (1973) proposed a degree of freedom correction for the M-estimate scale parameter.  The correction is given by \eqn{K = 1 + [p*(1-h_c)/n*h_c]} where \eqn{h_c} is the proportion of standardized residuals in absolute value less than the parameter \code{hparm}.  This correction \eqn{K} is used as a multiplicative factor to tauhat.  The default value of hparm is set at 2.

The usual degrees of freedom correction, \eqn{\sqrt{n/(n-p)}}, is also used as a multiplicative factor to tauhat.

}
\value{
Length one numeric object.
}
\references{ 
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Huber, P.J. (1973), Robust regression: Asymptotics, conjectures and Monte Carlo, \emph{Annals of Statistics}, 1, 799--821.

Koul, H.L., Sievers, G.L., and McKean, J.W. (1987), An estimator of the scale parameter for the rank analysis of linear models under general score functions, \emph{Scandinavian Journal of Statistics}, 14, 131--141.

McKean, J. W. and Sheather, S. J. (1991), Small Sample Properties of Robust Analyses of Linear Models Based on R-Estimates: A Survey, in \emph{Directions in Robust Statistics and Diagnostics, Part II}, Editors: W.\ Stahel and S.\ Weisberg, Springer-Verlag: New York, 1--19.
}
\author{ Joseph McKean, John Kloke}
\seealso{ \code{\link{rfit}} }

\examples{
#  For a standard normal distribution the parameter tau has the value 1.023327 (sqrt(pi/3)).
set.seed(283643659)
n <- 12; p <- 6; y <- rnorm(n); x <- matrix(rnorm(n*p),ncol=p)
tau1 <- rfit(y~x)$tauhat; tau2 <- rfit(y~x,delta=0.95)$tauhat
c(tau1,tau2) # 0.5516708 1.0138415
n <- 120; p <- 6; y <- rnorm(n); x <- matrix(rnorm(n*p),ncol=p)
tau3 <- rfit(y~x)$tauhat; tau4 <- rfit(y~x,delta=0.95)$tauhat
c(tau3,tau4) # 1.053974 1.041783
}

