toxbdry <- function(pLo, pHi, n, cP0=0.1, cP1=0.9, ngrid=6, niter=5) {
  nlook <- length(n)
  ptox <- seq(pLo, pHi, length=ngrid)
  alpha0 <- cP0/nlook
  r0 <- qbinom(1-alpha0, n, pLo)
#  print(r0)
  pstop0 <- bdrycross.prob(n, r0, ptox)
  alpha1 <- cP0
  r1 <- qbinom(1-alpha1, n, pLo)
#  print(r1)
  pstop1 <- bdrycross.prob(n, r1, ptox)
  while(pstop1[1,2] < cP0) {
    alpha1 <- 1.4*alpha1
    r1 <- qbinom(1 - alpha1, n, pLo)
    pstop1 <- bdrycross.prob(n, r1, ptox)
  }
  iter <- 0
  while (max(r0 - r1) > 1 | iter <= niter) {
    iter <- iter + 1
    alpha <- alpha0 + (alpha1-alpha0)*(cP0-pstop0[1,2])/(pstop1[1,2]-pstop0[1,2])
    r <- qbinom(1-alpha, n, pLo)
#    print(alpha)
#    print(r)
    pstop <- bdrycross.prob(n, r, ptox)
    if (pstop[1,2] < cP0) {
      alpha0 <- alpha
      r0 <- r
      pstop0 <- pstop
    } else {
      alpha1 <- alpha
      r1 <- r
      pstop1 <- pstop
    }
  }
# combine the operating characteristics of low and high boundaries  
  bdry.oc <- cbind(pstop1, pstop0[,2:4])
  colnames(bdry.oc)[2:7] <- c("pcross.lo", "pstop.lo", "ess.lo", "pcross.hi", "pstop.hi", "ess.hi")
  if(pstop1[ngrid,2] < cP1 & pstop0[ngrid,2] < cP1) warning(paste("Max sample size", n[nlook], "may be small for the specified stopping probabilities\n"))
  out <- list("looks"=n, "lo.bdry"=r1, "hi.bdry"=r0, "bdry.oc"=bdry.oc, "bdry.alpha"=c(alpha1, alpha0))
  class(out) <- "toxbdry"
  out
}

print.toxbdry <- function(x, ...) {
# convert the boundary to human readable form
  n <- x$looks
  if (max(diff(n)) == 1) {
    ii <- c(which(diff(x$lo.bdry) > 0), length(n))
    bdry.lo <- paste(x$lo.bdry[ii], n[ii], sep="/")
    ii <- c(which(diff(x$hi.bdry) > 0), length(n))
    bdry.hi <- paste(x$hi.bdry[ii], n[ii], sep="/")
  } else {
    bdry.lo <- paste(x$lo.bdry, n, sep="/")
    bdry.hi <- paste(x$hi.bdry, n, sep="/")
  }
  cat("\n Toxicity boundary based on repeated significance testing \n\n")
  cat(" ******************************************************************\n")
  cat(" * Stop if the number of toxicities exceeds (i.e. >) the boundary *\n")
  cat(" ******************************************************************\n")
  cat("\n  Low boundary:",  bdry.lo, "\n", sep="   ")
  cat(" High boundary:",  bdry.hi, "\n", sep="   ")
  cat("\n Operating Characteristics: \n\n")
  bdry.oc <- round(x$bdry.oc, digits=3)
  bdry.oc[,c(4,7)] <- round(bdry.oc[,c(4,7)], digits=1)
  print(bdry.oc)
}

bdrycross.prob <- function(n, r, ptox) {
  nlook <- length(n)
  np <- length(ptox)
  dn <- diff(c(0,n))
  pcur <- matrix(dbinom(rep(0:r[1],np), dn[1], rep(ptox, rep(r[1]+1, np))), r[1]+1, np)
  pstop <- 1-pbinom(r[1], dn[1], ptox)
  ess <- n[1]*pstop
  for(i in 2:nlook) {
    ll <- r[i-1] + 1
    pstop0 <- rep(0,np)
    for(j in (r[i]-r[i-1]):min(dn[i],r[i])) {
      pstop0 <- pstop0 + pcur[ll,]*(1-pbinom(j, dn[i], ptox))
      ll <- ll - 1
    }
    pstop <- pstop + pstop0
    ess <- ess + n[i]*pstop0
    pnext <- matrix(0, r[i]+1, np)
    for(j in 0:r[i]) {
      for(k in max(0, j-dn[i]):min(j,r[i-1])) {
        pnext[j+1,] <- pnext[j+1,] + pcur[k+1,]*dbinom(j-k, dn[i], ptox)
      }
    }
    pcur <- pnext
  }
  ess <- ess + n[nlook]*(1-pstop)
  pcross <- pstop
  pstop <- pstop - pstop0
  cbind(ptox, pcross, pstop, ess)
}