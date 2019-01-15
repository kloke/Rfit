wald.test.overall <- function(fit) {

  p <- length(coef(fit))
  betahat <- coef(fit)[-1]
  v <- vcov(fit)[2:p,2:p]

  df1 <- p - 1
  teststat <- t(betahat)%*%solve(v)%*%betahat/df1
  df2 <- length(fit$y) - fit$qrx1$rank

  pval <- pf(teststat,df1,df2,lower.tail=FALSE)

  ans <- list(F=teststat,p.value=pval,df1=df1,df2=df2)
  ans 

}
