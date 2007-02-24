coxphQuantile <- function(phfit, xrange, p=0.5, whichx=1, otherx=NULL,
                           ngrid=256, ...) {
  if (class(phfit) != "coxph") stop("phfit shoud be coxph class object")
  cvtmean <- phfit$means
  loghr <- phfit$coef
  S0 <- survfit(phfit)
  stime <- S0$time
  ssurv <- S0$surv
  if (!missing(otherx)) {
    ssurv <- ssurv^(exp(sum(loghr[-whichx]*(otherx-cvtmean[-whichx]))))
  }
  xseq <- seq(xrange[1], xrange[2], length.out=ngrid)
  ii <- sapply((xseq - cvtmean[whichx])*loghr[whichx], function(xi, ssurv, p) {
    which(ssurv^(exp(xi)) < p)[1]
  }, ssurv, p)
  lines(xseq, stime[ii], ...)
}
