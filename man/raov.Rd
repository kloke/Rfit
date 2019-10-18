\name{raov}
\alias{raov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
R ANOVA 
}
\description{
Returns full model fit and robust ANOVA table for all main effects and interactions. 
}
\usage{
raov(f, data = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
an object of class formula
}
  \item{data}{
an optional data frame
}
  \item{\dots}{
additional arguments
}
}
\details{
Based on reduction in dispersion tests for testing main effects and interaction.
Uses an algorithm described in Hocking (1985).
}
\value{
\item{table}{Description of 'comp1'}
\item{fit}{full model fit returned from rfit}
\item{residuals}{the residuals, i.e. y-yhat}
\item{fitted.values}{ yhat = x betahat}
\item{call}{Call to the function}
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Hocking, R. R. (1985), \emph{The Analysis of Linear Models}, Monterey, California: Brooks/Cole.
}
\author{
Joseph McKean, John Kloke
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rfit}}, \code{\link{oneway.rfit}}
}
\examples{
raov(logSurv~Poison+Treatment,data=BoxCox)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ nonparametric }
\keyword{ robust }
