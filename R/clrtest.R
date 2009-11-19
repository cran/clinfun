uclrtest <- function(tim, sts, grp, cvt, bb=NULL) {
  ord <- order(tim, -sts)
  n <- length(tim)
  ng <- length(unique(grp))
  if (is.matrix(cvt)) {
    p <- ncol(cvt)
  } else {
    p <- 1
  }
  osts <- sts[ord]
  ogrp <- as.numeric(as.factor(grp[ord]))
  ocov <- cvt[ord]
  zz <- .Fortran("uclrst",
                 as.integer(n),
                 as.integer(ng),
                 as.integer(p),
                 as.double(osts),
                 as.integer(ogrp),
                 as.double(ocov),
                 a0=double(n),
                 a1=double(n*ng),
                 xi=double(p),
                 xj=double(p),
                 Vii=double(ng),
                 Vij=double(ng),
                 Vji=double(ng),
                 Vjj=double(ng),
                 Vidot=double(ng*n),
                 Vdotj=double(ng*n),
                 Vijij=double(ng*ng),
                 igrp=double(ng),
                 jgrp=double(ng),
                 lrmn=double(ng),
                 lrvar=double(ng*ng),
                 as.double(bb))
  print(c(zz$igrp, zz$jgrp, zz$igrp-zz$jgrp))
  ustat <- zz$lrmn
  uvar <- matrix(zz$lrvar,ng,ng)
  chi <- sum(solve(uvar[-1,-1], ustat[-1]) * ustat[-1])
  list(statistic=ustat, var=uvar, chisq=chi)
}
