\name{disp}
\alias{disp}
\title{ Jaeckel's Dispersion Function }
\description{
Returns the value of Jaeckel's dispersion function for given values of the regression coefficents.
}
\usage{
disp(beta, x, y, scores)
}
\arguments{
  \item{beta}{ p by 1 vector of regression coefficents }
  \item{x}{ n by p design matrix }
  \item{y}{ n by 1 response vector }
  \item{scores}{ an object of class scores }
}
\details{
Returns the value of Jaeckel's disperion function evaluated 
at the value of the parameters in the function call.
That is,
sum_{i=1}^n a(R(e_i)) * e_i where 
R denotes rank
and a(1) <= a(2) <= ... <= a(n) are the scores.
The residuals (e_i i=1,...n) are calculated y - x beta.
}
\references{ 
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Jaeckel, L. A. (1972). Estimating regression coefficients by minimizing the dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449 - 1458.
}
\author{John Kloke, Joseph McKean}


\seealso{
\code{\link{rfit}}
\code{\link{drop.test}}
\code{\link{summary.rfit}}
}

