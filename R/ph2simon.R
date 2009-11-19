ph2simon <- function(pu, pa, ep1, ep2, nmax = 100) {
  if(nmax > 1000) stop("nmax cannot exceed 1000")
  nmax1 <- nmax + 1
  m <- (nmax * (nmax + 3))/2
  n <- rep(1:nmax, (1:nmax) + 1)
  x <- n * 0
  j <- 1
  for(i in 1:nmax) {
    x[j + (0:i)] <- 0:i
    j <- j + i + 1
  }
  p0 <- dbinom(x, n, pu)
  p1 <- dbinom(x, n, pa)
  cp0 <- 1 - pbinom(x - 1, n, pu)
  cp1 <- 1 - pbinom(x - 1, n, pa)
  zz <- .Fortran("f2bdry",
           as.integer(m),
           as.integer(nmax),
           ep1,
           ep2,
           p0,
           p1,
           cp0,
           cp1,
           bdry = as.integer(rep(0, nmax*4)),
           peten = rep(0, nmax*2),
           as.integer(nmax1),
           rep(0,nmax1),
           rep(0,nmax1),
           PACKAGE="clinfun")
  ph2 <- list()
  ph2out <- cbind(matrix(zz$bdry,nmax,4),matrix(zz$peten,nmax,2))
  ph2out <- ph2out[ph2out[,5]!=0,]
  if(is.null(dim(ph2out))) {
    errmesg <- paste("  No feasible solution found. \n\tIncrease maximum sample size.  Current nmax value = ",nmax,".",sep="")
    stop(message=errmesg)
  }
  dimnames(ph2out) <- list(NULL,c("r1", "n1", "r", "n", "EN(p0)", "PET(p0)"))
  ph2$pu <- pu
  ph2$pa <- pa
  ph2$alpha <- ep1
  ph2$beta <- ep2
  ph2$out <- ph2out
  ph2$nmax <- nmax
  class(ph2) <- "ph2simon"
  ph2
}

print.ph2simon <- function(x, ...) {
  xout <- x$out
  nmax <- x$nmax
  n <- nrow(xout)
  nopt <- ((1:n)[xout[,5]==min(xout[,5])])[1]
  xopt <- xout[c(nopt,1),]
  dimnames(xopt)[[1]] <- c("Optimal","Minimax")
  cat("\n Simon 2-stage Phase II design \n\n")
  cat("Unacceptable response rate: ",x$pu,"\n")
  cat("Desirable response rate: ",x$pa,"\n")
  cat("Error rates: alpha = ",x$alpha,"; beta = ",x$beta,"\n\n")
  print(xopt, digits = 4, ...)
  cat("\n")
  if(xopt[1,4]>nmax-10) warning(paste("  Optimal sample size too close to nmax. \n  Try increasing nmax (current value = ",nmax,")\n",sep=""))
}

plot.ph2simon <- function(x, ...) {
  xout <- x$out
  n <- nrow(xout)
  nopt <- ((1:n)[xout[,5]==min(xout[,5])])[1]
  nopt1 <- min(nopt+5,n)
  plot(xout[1:nopt1,4],xout[1:nopt1,5],type="l",xlab="Maximum number of patients",ylab="Expected trial size", ...)
  points(xout[1,4],xout[1,5],pch="M")
  points(xout[nopt,4],xout[nopt,5],pch="O")
}

oc.twostage.bdry <- function(pu, pa, r1, n1, r, n){
  pet <- err0 <- pbinom(r1,n1,pu)
  ess <- n1 + (n-n1)*(1-pet)
  err1 <- pbinom(r1,n1,pa)
  for(i in (r1+1):r) {
    err0 <- err0 + dbinom(i,n1,pu)*pbinom(r-i,n-n1,pu)
    err1 <- err1 + dbinom(i,n1,pa)*pbinom(r-i,n-n1,pa)
  }
  out <- c(1-err0, 1-err1, pet, ess)
  names(out) <- c("P(reject H0 | p0)","P(reject H0 | p1)","PET(p0)","EN(p0)")
  out
}
