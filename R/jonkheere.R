# The exact/permutation version of Jonckheere-Terpstra test
# asymptotic version is equivalent to cor.test for Kendall's tau

# this function computes the pdf of the statistic for small samples
djonckheere <- function(gsize) {
  ng <- length(gsize)
  cgsize <- rev(cumsum(rev(gsize)))
  mxsum <- sum(gsize[-ng]*cgsize[-1]) + 1
  jrsum <- rep(0, mxsum)
  jrsum[1] <- 1
  zz <- .Fortran("djonck",
                 as.integer(mxsum),
                 jrsum=as.double(jrsum),
                 as.integer(ng),
                 as.integer(cgsize))
  zz$jrsum/sum(zz$jrsum)
}

# now do the test as in wilcox.test
jonckheere.test <- function(x, g, alternative = c("two.sided", "increasing", "decreasing")) {
  if(!is.numeric(x)) stop("data values should be numeric")
  if(!is.numeric(g) & !is.ordered(g)) stop("group should be numeric or ordered factor")
  alternative <- match.arg(alternative)
  METHOD <- "Jonckheere-Terpstra test"
  n <- length(x)
  if(length(g) != n) stop("lengths of data values and group don't match")
  TIES <- length(unique(x)) != n
  gsize <- table(g)
  ng <- length(gsize)
  cgsize <- c(0,cumsum(gsize))
  x <- x[order(g)]
  jtrsum <- jtmean <- jtvar <- 0
  for(i in 1:(ng-1)) {
    na <- gsize[i]
    nb <- n-cgsize[i+1]
    jtrsum <- jtrsum + sum(rank(x[(cgsize[i]+1):n])[1:na]) - na*(na+1)/2
    jtmean <- jtmean + na*nb/2
    jtvar <- jtvar + na*nb*(na+nb+1)/12
  }
# jtrsum will be small if data are increasing and large if decreasing
  STATISTIC <- jtrsum
  names(STATISTIC) <- "JT"
  if (TIES | (n > 50)) {
    zstat <- (STATISTIC-jtmean)/sqrt(jtvar)
    PVAL <- pnorm(zstat)
    PVAL <- switch(alternative,
                   "two.sided" = 2*min(PVAL, 1-PVAL),
                   "increasing" = PVAL,
                   "decreasing" = 1-PVAL)
  } else {
    iPVAL <- sum(djonckheere(gsize)[1:(jtrsum+1)])
    dPVAL <- 1-sum(djonckheere(gsize)[1:(jtrsum)])
    PVAL <- switch(alternative,
                   "two.sided" = 2*min(iPVAL, dPVAL),
                   "increasing" = iPVAL,
                   "decreasing" = dPVAL)
  }
  RVAL <- list(statistic = STATISTIC,
               p.value = as.numeric(PVAL),
               alternative = alternative,
               method = METHOD)
  class(RVAL) <- "htest"
  RVAL
}
