print.summary.rfit <- function (x, digits = max(5, .Options$digits - 2), ...) {
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    est <- x$coef
    printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
    if( x$overall.test == 'drop') { 
    cat("\nMultiple R-squared (Robust):", x$R2, "\n")
    cat("Reduction in Dispersion Test:", round(x$dropstat, digits = digits), 
        "p-value:", round(x$droppval, digits = digits), "\n")
	} 
	if( x$overall.test == 'wald') {
    cat("Overall Wald Test:", round(x$waldstat, digits=digits),
        "p-value:", round(x$waldpval, digits = digits), "\n")
	} 
    cat("\n")
}
