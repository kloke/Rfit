\name{BoxCox}
\alias{BoxCox}
\docType{data}
\title{
Box and Cox (1964) data.
}
\description{
The data are  the results of a 3 * 4 two-way design, where
forty-eight animals were exposed to three different poisons and four different treatments.
The design is balanced with four replications per cell.
The response was the log survival time of the animal.
}
\usage{data(BoxCox)}
\format{
  A data frame with 48 observations on the following 3 variables.
  \describe{
    \item{\code{logSurv}}{log Survival Time}
    \item{\code{Poison}}{a factor indicating poison level}
    \item{\code{Treatment}}{a factor indicating treatment level}
  }
}
\source{
Box, G.E.P. and Cox, D.R. (1964), An analysis of transformations, \emph{ Journal of the Royal Statistical Society, Series B, Methodological}, 26, 211-252.
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\examples{
data(BoxCox)
with(BoxCox,interaction.plot(Treatment,Poison,logSurv,median))
raov(logSurv~Poison+Treatment,data=BoxCox)
}
\keyword{datasets}
