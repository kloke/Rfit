\name{CardioRiskFactors}
\alias{CardioRiskFactors}
\docType{data}
\title{
Cardiovascular risk factors
}
\description{
Data from a study to investigate assocation between uric acid and various cardiovascular risk factors in developing countries (Heritier et. al. 2009).
There are 474 men and 524 women aged 25-64.
}
\usage{data(CardioRiskFactors)}
\format{
  A data frame with 998 observations on the following 14 variables.
  \describe{
    \item{\code{age}}{Age of subject}
    \item{\code{bmi}}{Body Mass Index}
    \item{\code{waisthip}}{waist/hip ratio(?)}
    \item{\code{smok}}{indicator for regular smoker}
    \item{\code{choles}}{total cholesterol}
    \item{\code{trig}}{triglycerides level in body fat}
    \item{\code{hdl}}{high-density lipoprotien(?)}
    \item{\code{ldl}}{low-density lipoprotein}
    \item{\code{sys}}{systolic blood pressure}
    \item{\code{dia}}{diastolic blood pressure(?)}
    \item{\code{Uric}}{serum uric}
    \item{\code{sex}}{indicator for male}
    \item{\code{alco}}{alcohol intake (mL/day)}
    \item{\code{apoa}}{apoprotein A}
  }
}
\details{
Data set and description taken from Heritier et. al. (2009) (c.f. Conen et. al. 2004).
Some not discussed (in Section 3.5) of the text and their persummed meaning is listed followed by (?).
}
\source{
Heritier, S., Cantoni, E., Copt, S., and Victoria-Feser, M. (2009), \emph{Robust Methods in Biostatistics}, New York: John Wiley \& Sons.

Conen, D., Wietlisbach, V., Bovet, P., Shamlaye, C., Riesen, W., Paccaud, F., and Burnier, M. (2004), Prevalence of hyperuricemia and relation of serum uric acid with cardiovascular risk factors in a developing country.
\emph{BMC Public Health},
\url{http://www.biomedcentral.com/1471-2458/4/9}.
}
%\references{
%}
\examples{
data(CardioRiskFactors)
fitF<-rfit(Uric~bmi+sys+choles+ldl+sex+smok+alco+apoa+trig+age,data=CardioRiskFactors)
fitR<-rfit(Uric~bmi+sys+choles+ldl+sex,data=CardioRiskFactors)
drop.test(fitF,fitR)
summary(fitR)
}
\keyword{datasets}
