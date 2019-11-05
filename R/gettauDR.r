gettauDR <- function(e, p=0, scores, B=1002, h=bw.nrd(e)) {

warning('using experimental function to estimate tau based on binned residuals.')

  get_breaks <- function(ehat,B) {
    eps <- 0.00001
    breaks <- quantile(ehat,seq(0,1,length=B))
    ind <- c(1,length(breaks))
    breaks[ind] <- breaks[ind] + eps*sd(breaks)*c(-1,1)
    unique(breaks)
  }

  getStepScoresDeriv <- function(scores,counts) {
    rank1 <- cumsum(counts)        #high rank
    rank0 <- rank1 - counts + 1    #low rank
    ave_rank <- (rank1 + rank0)/2
    getScoresDeriv(scores,ave_rank/(sum(counts)+1))
  }

# Huber's DF correction
dftauhat <- function(tauhat,ehat,param,p){

        epshc <- 0.000001
        n <- length(ehat)
  hubcor <- mean( abs(ehat/tauhat) < param )
        if(hubcor < epshc){hubcor <- epshc}
        corr <- sqrt(n/(n-p))*(1 + ((p/n)*((1-hubcor)/hubcor)))
        tauhat <- tauhat*corr
        tauhat
}



  ngb <- hist(e,br=get_breaks(e,B),plot=FALSE)

  scrDs <- getStepScoresDeriv(scores=scores,counts=ngb$counts)

  D <- dnorm( outer(ngb$mids,ngb$mids,'-')/h,0.0,1.0,0 )/h
  C <- outer( ngb$counts, ngb$counts,'*')
  est <- drop(crossprod( apply(C*D,1,sum), scrDs ))

  n <- length(e)

  est <- est/n/n
  est <- 1/est

huber <- TRUE
if(huber) est <- dftauhat(est,e,param=2,p=p)

est


}
