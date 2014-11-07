\name{oneway.rfit}
\alias{oneway.rfit}
%\alias{summary.oneway.rfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rank-based Oneway Analysis of Variance }
\description{
Carries out a robust analysis of variance for a one factor design.
Analysis is based on the R estimates.
}
\usage{
oneway.rfit(y, g, scores = Rfit::wscores, p.adjust = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ n by 1 response vector }
  \item{g}{ n by 1 vector representing group membership }
  \item{scores}{ an object of class 'scores'}
  \item{p.adjust}{ adjustment to the p-values, argument passed to p.adjust }
}
\details{
Carries out a robust one-way analysis of variance based on full model r fit.
}
\value{
\item{ fit }{ full model fit from rfit }
\item{ est }{ Estimates }
\item{ se }{ Standard Errors } 
\item{ I }{ First Index }
\item{ J }{ Second Index }
\item{p.value}{ p-values }
\item{y}{ response vector }
\item{g}{ vector denoting group membership }
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\author{ Joseph McKean, John Kloke}

\seealso{\link{rfit} }
\examples{
	data(quail)
	oneway.rfit(quail$ldl,quail$treat)
 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ anova }
\keyword{ nonparametric }
\keyword{ robust }
