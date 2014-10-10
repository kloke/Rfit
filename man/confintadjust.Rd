\name{confintadjust}
\alias{confintadjust}
\alias{confintadjust.methods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence interval adjustment methods
}
\description{
Returns the critical value to be used in calculating adjusted confidence intervals.
Currently provides methods for Boneferroni and Tukey for confidence interval adjustment methods as well as no adjustment.
}
\usage{
confintadjust(n, k, alpha = 0.05, method = confintadjust.methods, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size
}
  \item{k}{
number of comparisons
}
  \item{alpha}{
overall (experimentwise) type I error rate
}
  \item{method}{
one of confintadjust.methods
}
  \item{\dots}{
Additonal arguments.  Currently not used.
}
}
\details{
Returns critial value based on one of the adjustment methods.
}
\value{
\item{cv}{critical value}
\item{method}{the method used}
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\author{
Joseph McKean, John Kloke
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(n,k,alpha=0.05,method=confintadjust.methods,...) {
	method<-match.arg(method)
	cv<-switch(method, bonferroni = qt(1-alpha/choose(k,2),n-k),
		tukey = qtukey(1-alpha,k,n-k)/sqrt(2),
		none = qt(1-alpha/2,n-k)
	)

	res<-list(cv=cv,method=method)
	res

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
