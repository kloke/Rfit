taustar <- function(e,p,conf=0.95) {

  n <- length(e)
  zc <- qnorm( (1+conf)/2)
  
  c1 <- floor(0.5*n - 0.5*sqrt(n)*zc - 0.5)

  if( c1 < 0 ) c1 <- 0

  ind <- c(c1+1,n-c1)
  z <- diff(sort(e,partial=ind)[ind])

  sqrt(n/(n-p-1)) * sqrt(n)*z/(2*zc)

}
