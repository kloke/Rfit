\name{grad}
\alias{grad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate the Gradiant of Jaekel's Dispersion Function }
\description{
 Calculate the Gradiant of Jaekel's Dispersion Function 
}
\usage{
grad(x, y, beta, scores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ n by p design matrix }
  \item{y}{ n by 1 response vector }
  \item{beta}{ p by 1 vector of regression coefficients} 
  \item{scores}{ an object of class scores }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The gradiant evaluated at beta.
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Jaeckel, L. A. (1972). Estimating regression coefficients by minimizing the dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449 - 1458.

Jureckova, J. (1971). Nonparametric estimate of regression coefficients. \emph{Annals of Mathematical Statistics}, 42, 1328 - 1338.
}
\author{ John Kloke \email{kloke@biostat.wisc.edu}}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## The function is currently defined as
function (x, y, beta, scores) 
{
    x <- as.matrix(x)
    e <- y - x \%*\% beta
    r <- rank(e, ties.method = "first")/(length(e) + 1)
    -t(x) \%*\% scores@phi(r)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
