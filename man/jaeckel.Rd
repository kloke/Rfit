\name{jaeckel}
\alias{jaeckel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to Minimize Jaeckel's Dispersion Function}
\description{
Uses the built-in function \code{optim} to minimize Jaeckel's dispersion function.  Alternates between \code{CG} and \code{BFGS} steps and initially takes a number of IRLS steps.
}
\usage{
jaeckel(x, y, beta0 = rq(y ~ x)$coef[2:(ncol(x) + 1)], scores = Rfit::wscores,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ n by p design matrix }
  \item{y}{ n by 1 response vector }
  \item{beta0}{ intial estimate }
  \item{scores}{ object of class 'scores' }
  \item{\dots}{addtional arguments to be passed to fitting routine }
}
\details{
This function is meant to mimic the minization algorithm implemented in RGLM (Kapenga, et. al. 1988). The main loop of the function alternates between CG and BFGS estimation methods.  To speed up convergence, it first takes a number of iterated reweighted least squares (IRLS) steps which generally gets close to the solution in a small number of steps.  
Using IRLS for rank regression was first considered by Cheng and Hettmansperger (1983).  See also Sievers and Abebe (2004).
}
\value{
Results of \code{optim} are returned.
}
\references{ 
Cheng, K. S. and Hettmansperger, T. P. (1983), Weighted Least-Squares Rank Regression, \emph{Communications in Statistics, Part A - Theory and Methods}, 12, 1069-1086.

Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Jaeckel, L. A. (1972), Estimating regression coefficients by minimizing the dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449 - 1458.

Kapenga, J. A., McKean, J. W., and Vidmar, T. J. (1988), \emph{RGLM: Users Manual}, Statist. Assoc. Short Course on Robust Statistical Procedures for the Analysis of Linear and Nonlinear Models, New Orleans.

Sievers, J. and Abebe, A. (2004), Rank Estimation of Regression Coefficients Using Iterated Reweighted Least Squares, \emph{Journal of Statistical Computation and Simulation}, 74, 821-831.

}
\author{ John Kloke \email{kloke@biostat.wisc.edu}}

\seealso{ \code{\link{optim}}, \code{\link{rfit}} }
\examples{
##  This is a internal function.  See rfit for user-level examples.

## The function is currently defined as
function (x, y, beta0 = rq(y ~ x - 1)$coef, scores = wscores, 
    maxiter = 100, irls0 = 10, BFGS0 = 20, stepCG = 5, stepBFGS = 2) 
{
    x <- x - outer(rep(1, nrow(x)), apply(x, 2, mean))
    beta0 <- irls(x, y, beta0, max.iter = irls0)
    if (BFGS0 < 1) 
        BFGS0 <- 1
    fit <- optim(beta0, disp, method = "BFGS", x = x, y = y, 
        scores = scores, control = list(maxit = BFGS0))
    iter <- 0
    while (fit$convergence && iter < maxiter) {
        iter <- iter + 1
        fit <- optim(fit$par, disp, method = "CG", x = x, y = y, 
            scores = scores, control = list(maxit = stepCG))
        fit <- optim(fit$par, disp, method = "BFGS", x = x, y = y, 
            scores = scores, control = list(maxit = stepBFGS))
    }
    optim(fit$par, disp, method = "BFGS", x = x, y = y, scores = scores)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
