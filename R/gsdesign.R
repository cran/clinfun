gsd.drift <- function(ifrac, sig.level = 0.05, pow = 0.8, delta.eb = 0.5, delta.fb = 0, alternative = c("two.sided", "one.sided"), tol=0.0001) {
  alternative <- match.arg(alternative)
  nlook <- length(ifrac)
# check that the last value of information fraction is 1
  if (ifrac[nlook] != 1) stop("last information fraction (ifrac) value should be 1")
# The statistics are standardized Brownian motion in [0,1] time interval
# The times are the proportion of information available.
# initialize the correlation matrix
  if (nlook == 1) {
    if (alternative=="one.sided") {
      efbdry <- qnorm(1-sig.level)
    } else {
      efbdry <- qnorm(1-sig.level/2)
    }
    drift0 <- efbdry + qnorm(pow)
  } else {
# check that information fraction is increasing
    if (!all(diff(ifrac) > 0)) stop("information fraction (ifrac) values should be increasing")
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
# iter <- 1
    if (alternative=="one.sided") {
      zlo <- qnorm(1-sig.level)
      zhi <- qnorm(1-sig.level/nlook)
      alo <- 1 - pmvnorm(upper=zlo/ebden, corr=corr)
      ahi <- 1 - pmvnorm(upper=zhi/ebden, corr=corr)
      zz <- zlo + (alo-sig.level)*(zhi-zlo)/(alo-ahi)
      aa <- 1 - pmvnorm(upper=zz/ebden, corr=corr)
# cat(paste("iteration =",iter,"; sig =",aa,"\n"))
      while(abs(aa-sig.level) > tol) {
# iter <- iter + 1
        if (aa>sig.level) {
          zlo <- zz
          alo <- aa
        } else {
          zhi <- zz
          ahi <- aa
        }
        zz <- zlo + (alo-sig.level)*(zhi-zlo)/(alo-ahi)
        aa <- 1 - pmvnorm(upper=zz/ebden, corr=corr)
# cat(paste("iteration =",iter,"; sig =",aa,"\n"))
      }
    } else {
      zlo <- qnorm(1-sig.level/2)
      zhi <- qnorm(1-sig.level/(2*nlook))
      alo <- 1 - pmvnorm(lower=-zlo/ebden, upper=zlo/ebden, corr=corr)
      ahi <- 1 - pmvnorm(lower=-zhi/ebden, upper=zhi/ebden, corr=corr)
      zz <- zlo + (alo-sig.level)*(zhi-zlo)/(alo-ahi)
      aa <- 1 - pmvnorm(lower=-zz/ebden, upper=zz/ebden, corr=corr)
# cat(paste("iteration =",iter,"; sig =",aa,"\n"))
      while(abs(aa-sig.level) > tol) {
# iter <- iter + 1
        if (aa>sig.level) {
          zlo <- zz
          alo <- aa
        } else {
          zhi <- zz
          ahi <- aa
        }
        zz <- zlo + (alo-sig.level)*(zhi-zlo)/(alo-ahi)
        aa <- 1 - pmvnorm(lower=-zz/ebden, upper=zz/ebden, corr=corr)
# cat(paste("iteration =",iter,"; sig =",aa,"\n"))
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
  }
  list("ifrac"=ifrac, "sig.level"=sig.level, "power"=pow, "alternative"=alternative, "delta.eb"=delta.eb, "efbdry"=efbdry, "drift"=drift0)
}


# drift parameter theta
# binomial: n (per arm) = theta^2 * 2*pbar*(1-pbar)/(pC -pE)^2
# normal:   n (per arm) = theta^2 * 2*sigma^2/(muC-muE)^2
# survival: d (total) = theta^2 * 4/(log(haz-ratio))^2
#           Convert number of events d to sample size n

gsdesign.binomial <- function(ifrac, pC, pE, sig.level=0.05, power=0.8,
                              delta.eb = 0.5, delta.fb = 0, alternative =
                              c("two.sided", "one.sided"), tol=0.0001) {
  drift.out <- gsd.drift(ifrac, sig.level, power, delta.eb, delta.fb, alternative, tol)
  pbar <- (pC+pE)/2
  n <- 2*pbar*(1-pbar)*(drift.out$drift/(pC-pE))^2
  out <- drift.out
  out$drift <- NULL
  out$pC <- pC
  out$pE <- pE
  out$outcome <- "binary"
  out$sample.size <- n
  class(out) <- "gsdesign"
  out
}

gsdesign.normal <- function(ifrac, delta, sd=1, sig.level=0.05, power=0.8,
                            delta.eb = 0.5, delta.fb = 0, alternative =
                            c("two.sided", "one.sided"), tol=0.0001) {
  drift.out <- gsd.drift(ifrac, sig.level, power, delta.eb, delta.fb, alternative, tol)
  n <- 2*(drift.out$drift*sd/delta)^2
  out <- drift.out
  out$drift <- NULL
  out$delta <- delta
  out$sd <- sd
  out$outcome <- "normal"
  out$sample.size <- n
  class(out) <- "gsdesign"
  out
}

gsdesign.survival <- function(ifrac, haz.ratio, sig.level = 0.05, power = 0.8,
                              delta.eb = 0.5, delta.fb = 0, alternative = 
                              c("two.sided", "one.sided"), tol=0.0001) {
  drift.out <- gsd.drift(ifrac, sig.level, power, delta.eb, delta.fb, alternative, tol)
  d <- 4*(drift.out$drift/log(haz.ratio))^2
  out <- drift.out
  out$drift <- NULL
  out$haz.ratio <- haz.ratio
  out$outcome <- "survival"
  out$num.events <- d
  class(out) <- "gsdesign"
  out
}

print.gsdesign <- function(x, ...) {
  if (class(x) != "gsdesign") stop("input shoud be a gsdesign class object")
  cat("\n Group sequential design for comparing", x$outcome, "data with ")
  switch(match(x$outcome, c("binary", "normal", "survival")),
         cat("rates  pC =", x$p1,", pE =", x$p2, "\n\n"),
         cat("delta =", x$delta, ", sd =", x$sd, "\n\n"),
         cat("hazard ratio =", x$haz.ratio, "\n\n"))
  switch(match(x$outcome, c("binary", "normal", "survival")),
         cat("  sample size (per arm) =", x$sample.size, "\n"),
         cat("  sample size (per arm) =", x$sample.size, "\n"),
         cat(" total number of events =", x$num.events, "\n"))
  cat("   information fraction =", format(round(x$ifrac, 3), digits=3), "\n")
  cat("          boundary type =", x$delta.eb, " (Pocock = 0; O'Brien-Fleming = 0.5) \n")
  cat("      efficacy boundary =", round(x$efbdry, 3), "\n")
  cat("              sig.level =", x$sig.level, "\n")
  cat("                  power =", x$power, "\n")
  cat("            alternative =", x$alternative, "\n\n")
  invisible(x)
}
