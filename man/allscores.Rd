\name{allscores}
\alias{wscores}
\alias{nscores}
\alias{bentscores1}
\alias{bentscores2}
\alias{bentscores3}
\alias{bentscores4}
\alias{logGFscores}
\alias{logrank.scores}
\docType{data}
\title{ All Scores }
\description{
An object of class scores which includes the score function and it's derivative for rank-based regression inference.
}
\usage{data(wscores)}
\format{
  The format is:
Formal class 'scores' [package ".GlobalEnv"] with 2 slots
  ..@ phi :function (u)  
  ..@ Dphi:function (u)  
}
\details{
Using Wilcoxon (linear) scores leads to inference which has ARE of 0.955 to least squares (ML) when the data are normal.
Wilcoxon scores are optimal when the underlying error distribution is logistic.
Normal scores are optimal when the data are normally distributed.
Log-rank scores are optimal when the data are from an exponential distribution, e.g. in a proportional hazards model.
Log-Generalized F scores can also be used in the analysis of survival data (see Hettmansperger and McKean p. 233).

bentscores1 are recommended for right-skewed distributions.
bentscores2 are recommended for light-tailed distributions.
bentscores3 are recommended for left-skewed distributions.
bentscores4 are recommended for heavy-tailed distributions.
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\examples{
u <- seq(0.01,0.99,by=0.01)
plot(u,getScores(wscores,u),type='l',main='Wilcoxon Scores')
plot(u,getScores(nscores,u),type='l',main='Normal Scores')

data(wscores)
x<-runif(50)
y<-rlogis(50)
rfit(y~x,scores=wscores)

x<-rnorm(50)
y<-rnorm(50)
rfit(y~x,scores=nscores)

}
\keyword{datasets}
