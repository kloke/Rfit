#' Confidence Interval Method for Objects of Class Rfit  
#' 
#' Implements a method for calculating confidence intervals for robust linear models of class Rfit.  
#' Based on Kloke and McKean (2012:58) and the general design of confint.lm.
#' 
#' @param object An object of class rfit
#' @param level The confidence level for the interval of interest. Defaults to 0.95
#' @author Gustavo A. Ballen \email{gaballench@gmail.com}
#' @return A matrix with the lower and upper limits of the confidence interval at the specified level for
#' the coefficients in the object of interest.
#' @seealso \link{rfit}
#' @references Kloke, J.D. & McKean, J.W. (2012) Rfit: Rank-based estimation for linear models. \emph{The R Journal 4, 57–64}
#' @keywords rfit confidence interval
#' @examples
#' 
#' 	data(telephone)
#'      model <- rfit(calls ~ year, data = telephone)
#'      confint(model)
#'  
#' @export confint.rfit
#require(complmrob)

# rfit method for confint
# function copied from stats
format.perc <- function (probs, digits) {
    paste(format(100 * probs, trim = TRUE, 
                 scientific = FALSE, digits = digits), "%")
}

confint.rfit <- function(object, level = 0.95) {
    alpha <- 1 - level
    t.value <- qt(p = (1 - alpha/2), df = (length(object$fitted.values) - length(object$coefficients) - 1))
    se.beta <- coef(summary(object))[, "Std. Error"]
    coefNames <- rownames(coef(summary(object)))
    percNames <- format.perc(c((1 - (1 - alpha/2)), (1 - alpha/2)), 3)
    output <- array(NA,
                    dim = c(length(object$coefficients), length(coefNames)),
                    dimnames = list(coefNames, percNames))
    output[, 1] <- coef(summary(object))[, "Estimate"] - t.value * se.beta
    output[, 2] <- coef(summary(object))[, "Estimate"] + t.value * se.beta
    return(output)
}
