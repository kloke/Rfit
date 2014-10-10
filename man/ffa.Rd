\name{ffa}
\alias{ffa}
\docType{data}
\title{
Free Fatty Acid Data
}
\description{
The response variable is level of free fatty acid in a sample of prepubescent boys. The explanatory variables are age (in months), weight (in lbs), and skin fold thickness.
}
\usage{data(ffa)}
\format{
A data frame with 41 rows and 4 columns.
\describe{
	\item{\code{age}}{ age in years }
	\item{\code{weight}}{ weight in lbs }
	\item{\code{skin}}{ skin fold thinkness }
	\item{\code{ffa}}{  free fatty acid }
}
}
\source{
Morrison, D.F. (1983), \emph{Applied Linear Statistical Models}, Englewood Cliffs, NJ:Prentice Hall.
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\examples{
data(ffa)
summary(rfit(ffa~age+weight+skin,data=ffa))  #using the default (Wilcoxon scores)
summary(rfit(ffa~age+weight+skin,data=ffa,scores=bentscores1))
}
\keyword{datasets}
