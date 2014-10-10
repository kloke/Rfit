\name{serumLH}
\alias{serumLH}
\docType{data}
\title{
Serum Level of luteinizing hormone (LH)
}
\description{
Hollander and Wolfe (1999) discuss a 2 by 5 factorial design for a study
to determine the effect of light on the release of luteinizing hormone (LH).
The factors in the design are: light regimes at two levels (constant light and 14 hours of light followed by 10 hours of darkness) and a luteinizing release factor (LRF) at 5 different dosage levels.
The response is the level of luteinizing hormone (LH), nanograms per ml of serum in
blood samples.
Sixty rats were put on test under these 10 treatment combinations, six rats per combination.

}
\usage{data(serumLH)}
\format{
  A data frame with 60 observations on the following 3 variables.
  \describe{
    \item{\code{serum}}{a numeric vector}
    \item{\code{light.regime}}{a factor with levels \code{Constant} \code{Intermittent}}
    \item{\code{LRF.dose}}{a factor with levels \code{0} \code{10} \code{1250} \code{250} \code{50}}
  }
}
\source{
Hollander, M. and Wolfe, D.A. (1999), \emph{Nonparametric Statistical Methods}, New York: Wiley.
}
\references{
Hollander, M. and Wolfe, D.A. (1999), \emph{Nonparametric Statistical Methods}, New York: Wiley.
}
\examples{
data(serumLH)
raov(serum~light.regime + LRF.dose + light.regime*LRF.dose, data = serumLH)
}
\keyword{datasets}
