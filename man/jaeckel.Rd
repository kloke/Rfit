\name{jaeckel}
\alias{jaeckel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to Minimize Jaeckel's Dispersion Function}
\description{
Uses the built-in function \code{optim} to minimize Jaeckel's dispersion function. 
}
\usage{
jaeckel(x, y, beta0 = lm(y ~ x)$coef[2:(ncol(x) + 1)], 
  scores = Rfit::wscores, control = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ n by p design matrix }
  \item{y}{ n by 1 response vector }
  \item{beta0}{ intial estimate }
  \item{scores}{ object of class 'scores' }
  \item{control}{ control passed to fitting routine }
  \item{\dots}{addtional arguments to be passed to fitting routine }
}
\details{
Function uses \code{optim} with method set to BFGS to minimize Jaeckel's dispersion function.
If control is not specified at the function call, the relative tolerance (reltol) is set to .Machine$double.eps^(3/4)
maximum number of iterations is set to 200.
See \code{optim}.
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
