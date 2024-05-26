\name{jaeckel}
\alias{jaeckel}
\title{ Function to Minimize Jaeckel's Dispersion Function}
\description{
Uses the built-in function \code{optim} to minimize Jaeckel's dispersion function with respect to beta. 
}
\usage{
jaeckel(x, y, beta0 = lm(y ~ x)$coef[2:(ncol(x) + 1)], 
  scores = Rfit::wscores, control = NULL,...)
}
\arguments{
  \item{x}{ n by p design matrix }
  \item{y}{ n by 1 response vector }
  \item{beta0}{ initial estimate of beta }
  \item{scores}{ object of class 'scores' }
  \item{control}{ control passed to fitting routine }
  \item{\dots}{addtional arguments to be passed to fitting routine }
}
\details{
Jaeckel's dispersion function (Jaeckel 1972) is a convex function which measures the distance between the observed responses \eqn{y} and the fitted values \eqn{x \beta}.
The dispersion function is a sum of the products of the residuals, \eqn{y - x \beta}, and the scored ranks of the residuals.  A rank-based fit minimizes the dispersion function; see McKean and Schrader (1980) and Kloke and McKean (2012) for discussion.
\code{jaeckel} uses \code{optim} with the method set to BFGS to minimize Jaeckel's dispersion function.
If control is not specified at the function call, the relative tolerance (reltol) is set to \code{.Machine$double.eps^(3/4)} and maximum number of iterations is set to 200.

\code{jaeckel} is intended to be an internal function.  See \code{rfit} for a general purpose function.
}
\value{
Results of \code{optim} are returned.
}
\references{ 

Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Jaeckel, L. A. (1972), Estimating regression coefficients by minimizing the dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449 - 1458.

Kapenga, J. A., McKean, J. W., and Vidmar, T. J. (1988), \emph{RGLM: Users Manual}, Statist. Assoc. Short Course on Robust Statistical Procedures for the Analysis of Linear and Nonlinear Models, New Orleans.

}
\author{ John Kloke }

\seealso{ \code{\link{optim}}, \code{\link{rfit}} }
\examples{
##  This is a internal function.  See rfit for user-level examples.
}
