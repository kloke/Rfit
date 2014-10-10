\name{taufuncs}
\alias{hstarreadyscr}
\alias{hstar}
\alias{looptau}
\alias{pairup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Internal Functions for Estimating tau }
\description{
These are internal functions used for calculating the scale parameter tau necessary for estimating the standard errors of coefficients for rank-regression.
}
\usage{
hstarreadyscr(ehat,asc,ascpr)
hstar(abdord, wtord, const, n, y) 
looptau(delta, abdord, wtord, const, n)
pairup(x,type="less") 
}
\references{ 
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Koul, H.L., Sievers, G.L., and McKean, J.W. (1987) An esimator of the scale parameter for the rank analysis of linear models under general score functions, \emph{Scandinavian Journal of Statistics}, 14, 131-141.
}
\author{ Joseph McKean, John Kloke}

\arguments{
  \item{ehat}{ Full model residals }
  \item{delta}{Window parameter (proportion) used in the Koul et al. estimator of tau.  Default value is 0.80.   If the ratio of sample size to number of regression parameters (n to p) is less than 5, larger values such as 0.90 to 0.95 are more approporiate.}
  \item{y}{Argument of function hstar}
  \item{abdord}{Ordered absolute differences of residuals}
  \item{wtord}{Standardized (by const) ordered absolute differences of residuals}
  \item{const}{Range of score function}
  \item{n}{Sample size}
  \item{x}{Argument for pairup}
  \item{type}{Argument for the function pairup}
  \item{asc}{ scores }
  \item{ascpr}{ derivative of the scores }
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

\seealso{ \code{\link{gettau}}, \code{\link{rfit}} }
