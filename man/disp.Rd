\name{disp}
\alias{disp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Jaeckel's Dispersion Function }
\description{
Returns the value of Jaeckel's dispersion function for given values of the regression coefficents.
}
\usage{
disp(beta, x, y, scores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{ p by 1 vector of regression coefficents }
  \item{x}{ n by p design matrix }
  \item{y}{ n by 1 response vector }
  \item{scores}{ an object of class scores }
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
\references{ 
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Jaeckel, L. A. (1972). Estimating regression coefficients by minimizing the dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449 - 1458.
}
\author{ John Kloke \email{kloke@biostat.wisc.edu}}

%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\seealso{
%~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{summary.rfit}}
}

\examples{

## The function is currently defined as
function (beta, x, y, scores) 
{
    x <- as.matrix(x)
    e <- y - x \%*\% beta
    r <- rank(e, ties.method = "first")/(length(e) + 1)
    scores@phi(r) \%*\% e
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
