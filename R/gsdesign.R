gsd.drift <- function(ifrac, sig.level = 0.05, pow = 0.8, delta.eb = 0.5, delta.fb = 0, alternative = c("two.sided", "one.sided"), tol=0.001) {
  alternative <- match.arg(alternative)
  nlook <- length(ifrac)
# The statistics are standardized Brownian motion in [0,1] time interval
# The times are the proportion of information available.
# initialize the correlation matrix
  corr <- diag(1,nlook)
  for(i in 2:nlook) {
    for(j in 1:(i-1)) {
      corr[i,j] <- min(ifrac[i],ifrac[j])/sqrt(ifrac[i]*ifrac[j])
    }
  }
# The boundary is of the for constant/(ifrac^delta.eb) where
# delta.eb =0 for Pocock and =0.5 for O'Brien-Fleming boundary
# Calculate the constant for a given significance level under H0
  ebden <- ifrac^delta.eb
#  iter <- 1
  if (alternative=="one.sided") {
    zlo <- qnorm(1-sig.level)
    zhi <- qnorm(1-sig.level/nlook)
    alo <- 1 - pmvnorm(upper=zlo/ebden, corr=corr)
    ahi <- 1 - pmvnorm(upper=zhi/ebden, corr=corr)
    zz <- zlo + (alo-sig.level)*(zhi-zlo)/(alo-ahi)
    aa <- 1 - pmvnorm(upper=zz/ebden, corr=corr)
#    cat(paste("iteration =",iter,"; sig =",aa,"\n"))
    while(abs(aa-sig.level) > tol) {
#      iter <- iter + 1
      if (aa>sig.level) {
        zlo <- zz
        alo <- aa
      } else {
        zhi <- zz
        ahi <- aa
      }
      zz <- zlo + (alo-sig.level)*(zhi-zlo)/(alo-ahi)
      aa <- 1 - pmvnorm(upper=zz/ebden, corr=corr)
#      cat(paste("iteration =",iter,"; sig =",aa,"\n"))
    }
  } else {
    zlo <- qnorm(1-sig.level/2)
    zhi <- qnorm(1-sig.level/(2*nlook))
    alo <- 1 - pmvnorm(lower=-zlo/ebden, upper=zlo/ebden, corr=corr)
    ahi <- 1 - pmvnorm(lower=-zhi/ebden, upper=zhi/ebden, corr=corr)
    zz <- zlo + (alo-sig.level)*(zhi-zlo)/(alo-ahi)
    aa <- 1 - pmvnorm(lower=-zz/ebden, upper=zz/ebden, corr=corr)
#    cat(paste("iteration =",iter,"; sig =",aa,"\n"))
    while(abs(aa-sig.level) > tol) {
#      iter <- iter + 1
      if (aa>sig.level) {
        zlo <- zz
        alo <- aa
      } else {
        zhi <- zz
        ahi <- aa
      }
      zz <- zlo + (alo-sig.level)*(zhi-zlo)/(alo-ahi)
      aa <- 1 - pmvnorm(lower=-zz/ebden, upper=zz/ebden, corr=corr)
#      cat(paste("iteration =",iter,"; sig =",aa,"\n"))
    }
  }
  efbdry <- zz/ebden
# Under the alternative the Brownian motion has drift mu.
# So the standardized statistics have means mu*sqrt(ifrac)
# Calculate drift for given power.
  if (alternative=="one.sided") {
    drift0 <- qnorm(pow) + qnorm(1-sig.level)
  } else {
    drift0 <- qnorm(pow) + qnorm(1-sig.level/2)
  }
  driftlo <- drift0*0.8
  powlo <- 1 - pmvnorm(upper=efbdry-driftlo*sqrt(ifrac), corr=corr)
  drifthi <- drift0*1.25
  powhi <- 1 - pmvnorm(upper=efbdry-drifthi*sqrt(ifrac), corr=corr)
  drift0 <- driftlo + (pow-powlo)*(drifthi-driftlo)/(powhi-powlo)
  pow0 <- 1 - pmvnorm(upper=efbdry-drift0*sqrt(ifrac), corr=corr)
#  iter <- 1
#  cat(paste("iteration =",iter,"; drift =",drift0,"; pow =",pow0,"\n"))
#  while (abs(pow0-pow)>tol & iter<50) {
  while (abs(pow0-pow)>tol) {
    if (pow0>pow) {
      drifthi <- drift0
      powhi <- pow0
    } else {
      driftlo <- drift0
      powlo <- pow0
    }
    drift0 <- driftlo + (pow-powlo)*(drifthi-driftlo)/(powhi-powlo)
    pow0 <- 1 - pmvnorm(upper=efbdry-drift0*sqrt(ifrac), corr=corr)
#    iter <- iter+1
#    cat(paste("iteration =",iter,"; drift =",drift0,"; pow =",pow0,"\n"))
  }
  attributes(drift0) <- NULL
  list("ifrac"=ifrac, "sig.level"=sig.level, "power"=pow, "alternative"=alternative, "delta.eb"=delta.eb, "efbdry"=efbdry, "drift"=drift0)
}


# drift parameter theta
# binomial: n (per arm) = theta^2 * 2*pbar*(1-pbar)/(pC -pE)^2
# normal:   n (per arm) = theta^2 * 2*sigma^2/(muC-muE)^2
# survival: d (total) = theta^2 * 4/(log(haz-ratio))^2
#           Convert number of events d to sample size n

gsdesign.binomial <- function(ifrac, pC, pE, sig.level=0.05, power=0.8,
                              delta.eb = 0.5, delta.fb = 0, alternative =
                              c("two.sided", "one.sided"), tol=0.001) {
  drift.out <- gsd.drift(ifrac, sig.level, power, delta.eb, delta.fb, alternative)
  pbar <- (pC+pE)/2
  n <- 2*pbar*(1-pbar)*(drift.out$drift/(pC-pE))^2
  out <- drift.out
  out$drift <- NULL
  out$sample.size <- n
  out
}

gsdesign.normal <- function(ifrac, delta, sd=1, sig.level=0.05, power=0.8,
                            delta.eb = 0.5, delta.fb = 0, alternative =
                            c("two.sided", "one.sided"), tol=0.001) {
  drift.out <- gsd.drift(ifrac, sig.level, power, delta.eb, delta.fb, alternative)
  n <- 2*(drift.out$drift*sd/delta)^2
  out <- drift.out
  out$drift <- NULL
  out$sample.size <- n
  out
}

gsdesign.survival <- function(ifrac, haz.ratio, sig.level = 0.05, power = 0.8,
                              delta.eb = 0.5, delta.fb = 0, alternative = 
                              c("two.sided", "one.sided"), tol=0.001) {
  drift.out <- gsd.drift(ifrac, sig.level, power, delta.eb, delta.fb, alternative)
  d <- 4*(drift.out$drift/log(haz.ratio))^2
  out <- drift.out
  out$drift <- NULL
  out$num.events <- d
  out
}
