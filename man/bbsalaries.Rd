\name{bbsalaries}
\alias{bbsalaries}
\docType{data}
\title{Baseball Salaries}
\description{
Salaries of 176 professional baseball players for the 1987 season.
}
\usage{data(bbsalaries)}
\format{
  A data frame with 176 observations on the following 8 variables.
  \describe{
    \item{\code{logYears}}{Log of the number of years experience}
    \item{\code{aveWins}}{Average wins per year}
    \item{\code{aveLosses}}{Average losses per year}
    \item{\code{era}}{Earned Run Average}
    \item{\code{aveGames}}{Average games pitched in per year}
    \item{\code{aveInnings}}{Average number of innings pitched per year}
    \item{\code{aveSaves}}{Average number of saves per year}
    \item{\code{logSalary}}{Log of the base salary in dollars}
  }
}
%\details{
  %~~ If necessary, more details than the __description__ above ~~
%}
\source{
http://lib.stat.cmu.edu/datasets/baseball.data
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\examples{
data(bbsalaries)
summary(rfit(logSalary~logYears+aveWins+aveLosses+era+aveGames+aveInnings+aveSaves,data=bbsalaries))
}
\keyword{datasets}
