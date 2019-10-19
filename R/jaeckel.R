#' Function to Minimize Jaeckel's Dispersion Function
#' 
#' Uses the built-in function \code{optim} to minimize Jaeckel's dispersion
#' function.  
#' 
#' @param x n by p design matrix
#' @param y n by 1 response vector
#' @param beta0 intial estimate
#' @param scores object of class 'scores'
#' @return Results of \code{optim} are returned.
#' @author John Kloke 
#' @seealso \code{\link{optim}}, \code{\link{rfit}}
#' 
#' @references
#' Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric
#' Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Jaeckel, L. A. (1972), Estimating regression coefficients by minimizing the
#' dispersion of residuals. \emph{Annals of Mathematical Statistics}, 43, 1449
#' - 1458.
#' 
#' Kapenga, J. A., McKean, J. W., and Vidmar, T. J. (1988), \emph{RGLM: Users
#' Manual}, Statist. Assoc. Short Course on Robust Statistical Procedures for
#' the Analysis of Linear and Nonlinear Models, New Orleans.
#' 
#' ##  This is a internal function.  See rfit for user-level examples.
#' 
#' @export jaeckel
jaeckel <- function (x, y, beta0 = lm(y ~ x)$coef[2:(ncol(x) + 1)], 
  scores = Rfit::wscores, control=NULL, ...) {
  
# set max iter and/or convergence tolerance
if(is.null(control)) {
  control <- list(maxit=200,reltol=.Machine$double.eps^(3/4))
} else {
  if(is.null(control$reltol)) {
    control$reltol<-.Machine$double.eps^(3/4)
  }
  if(is.null(control$maxit)) {
     control$maxit<-200
  }
}

  scrs <- getScores(scores,seq_len(length(y))/(length(y)+1))

  j.grad <- function (x, y, beta, scrs,...) {
    x <- as.matrix(x)
    e <- y - x %*% beta
#    r <- rank(e, ties.method = "first")/(length(e) + 1)
#    crossprod(x,-1*getScores(scores,r) )
    crossprod(x[order(e),],-1*scrs )
  }


  j.disp <- function (beta, x, y,scrs,...) {
    e <- y - x %*% beta
    drop(crossprod(e[order(e)],scrs))
  }

  sd.y <- sd(y)
  ystar <- y/sd.y

  fit0 <- optim(beta0/sd.y, j.disp, method = "BFGS", x = x, y = ystar, scrs=scrs,
    gr=j.grad,control=control,...)

  optim(fit0$par*sd.y, j.disp, method = "BFGS", x = x, y = y, scrs=scrs,
    gr=j.grad,control=control,...)

}

