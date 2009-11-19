toxbdry2 <- function(pLo, pHi, n, cP0=0.1, cP1=0.9, ngrid=6, niter=5) {
# decide whether to satisfy null or the alternative
  nlook <- length(n)
  ptox <- rev(seq(pLo, pHi, length=ngrid))
# error threshold is cP0 for null and cP1 for alt
  alpha0 <- cP0/nlook
  r0 <- qbinom(1-alpha0, n, pLo)
#  print(r0)
  pstop0 <- bdrycross.prob(n, r0, ptox)
  alpha1 <- cP0
  r1 <- qbinom(1-alpha1, n, pLo)
#  print(r1)
  pstop1 <- bdrycross.prob(n, r1, ptox)
# depending on null or alt prob of interest is 1st or ngrid element
  while(pstop1[1,2] < cP1) {
    alpha0 <- alpha1
    pstop0 <- pstop1
    alpha1 <- 1.4*alpha1
    r1 <- qbinom(1 - alpha1, n, pLo)
    pstop1 <- bdrycross.prob(n, r1, ptox)
  }
  iter <- 0
  while (max(r0 - r1) > 1 | iter <= niter) {
    iter <- iter + 1
    if (alpha1 - alpha0 > 0.01) {
      alpha <- (alpha1 + alpha0)/2
    } else {
      alpha <- alpha0 + (alpha1-alpha0)*(cP1-pstop0[1,2])/(pstop1[1,2]-pstop0[1,2])
    }
    print(iter)
    r <- qbinom(1-alpha, n, pLo)
#    print(alpha)
#    print(r)
    pstop <- bdrycross.prob(n, r, ptox)
    if (pstop[1,2] < cP1) {
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
