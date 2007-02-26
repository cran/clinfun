# permlogrank computes the p-values based on permutation reference distribution
# 
# adapted from the original survdiff.s
permlogrank <- function(formula, data, subset, na.action, rho=0, nperm=5000) {
  call <- match.call()
  m <- match.call(expand=FALSE)
  m$rho <- NULL
  m$nperm <- NULL

  Terms <- if(missing(data)) terms(formula, 'strata')
  else              terms(formula, 'strata', data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  y <- model.extract(m, response)
  if (!inherits(y, "Surv")) stop("Response must be a survival object")
  if (attr(y, 'type') != 'right') stop("Right censored data only")
  ny <- ncol(y)
  n <- nrow(y)

  offset<- attr(Terms, "offset")
  if (!is.null(offset)) {
	#one sample test
    stop(" permutation logrank is for k-sample test only")
  }

  else { #k sample test
    strats <- attr(Terms, "specials")$strata
    if (length(strats)) {
      temp <- untangle.specials(Terms, 'strata', 1)
      dropx <- temp$terms
      if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
      else strata.keep <- strata(m[,temp$vars], shortlabel=TRUE)
    }
    else strata.keep <- rep(1,nrow(m))

    #Now create the group variable
    if (length(strats)) ll <- attr(Terms[-dropx], 'term.labels')
    else                ll <- attr(Terms, 'term.labels')
    if (length(ll) == 0) stop("No groups to test")
    else groups <- strata(m[ll])

    # re-order data
    ngroup <- length(unique(groups))
    x <- as.numeric(groups)
#    if (is.category(groups)) x <- as.numeric(groups)
#    else x <- match(groups, unique(groups))
    if (length(x) != n) stop("Data length mismatch")

    strat <- as.numeric(as.factor(strata.keep))
    nstrat <- length(unique(strat))
    if (length(strat) !=n) stop("Data length mismatch")

    ord <- order(strat, y[,1], -y[,2])
    y <- y[ord,]
    x <- x[ord]
    strat <- strat[ord]
    strat2 <- c(1*(diff(strat)!=0), 1)
    nstrat2 <- c(0,which(strat2==1))

    fit <- logrank.fit(n, ngroup, nstrat, rho, y, x, strat2)
    ochi <- fit$chi
    pchi <- rep(0,nperm)
    uu <- rep(0, n)
    for(i in 1:nperm) {
      uu <- runif(n)
      ord <- stratperm1(n, 1:n, nstrat, nstrat2, uu)
      pchi[i] <- logrank.fit(n, ngroup, nstrat, rho, y, x[ord], strat2)$chi
    }
    perm.p <- sum(ochi <= pchi)/nperm
    rval <-list(n= table(groups), obs = fit$observed,
                exp = fit$expected, var=fit$var,  chisq=fit$chi,
                perm.p = perm.p)
    if (length(strats)) rval$strata <- table(strata.keep)
  }

  na.action <- attr(m, "na.action")
  if (length(na.action)) rval$na.action <- na.action
  rval$call <- call
  class(rval) <- 'permlogrank'
  rval
}

# adapted from the orginal survdiff.fit.s
logrank.fit <- function(n, ngroup, nstrat, rho, y, x, strat2) {
  xx <- .C("survdiff2", as.integer(n),
           as.integer(ngroup),
           as.integer(nstrat),
           as.double(rho),
           as.double(y[,1]),
           as.integer(y[,2]),
           as.integer(x),
           as.integer(strat2),
           observed = double(ngroup*nstrat),
           expected = double(ngroup*nstrat),
           var.e    = double(ngroup * ngroup),
           double(ngroup), double(n),
           PACKAGE="survival")

  if (nstrat==1)  fit <- list(expected = xx$expected,
                              observed = xx$observed,
                              var      = matrix(xx$var, ngroup, ngroup))
  else            fit <- list(expected = matrix(xx$expected, ngroup),
			      observed = matrix(xx$observed, ngroup),
			      var      = matrix(xx$var, ngroup, ngroup))
  if (is.matrix(fit$observed)){
    otmp <- apply(fit$observed,1,sum)
    etmp <- apply(fit$expected,1,sum)
  }
  else {
    otmp <- fit$observed
    etmp <- fit$expected
  }
  df   <- (etmp >0)                #remove groups with exp=0
  if (sum(df) <2) chi <- 0         # No test, actually
  else {
    temp2 <- ((otmp - etmp)[df])[-1]
    vv <- (fit$var[df,df])[-1,-1, drop=FALSE]
    chi <- sum(solve(vv, temp2) * temp2)
  }
  fit$chi <- chi
  fit
}

# adapted from the original print.survdiff.s
print.permlogrank <- function(x, digits = max(options()$digits - 4, 3), ...) {

  saveopt <-options(digits=digits)
  on.exit(options(saveopt))

  if (!inherits(x, 'permlogrank'))
    stop("Object is not the result of permlogrank")
  if (!is.null(cl<- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }

  omit <- x$na.action
  if (length(omit)) cat("n=", sum(x$n), ", ", naprint(omit),
					  ".\n\n", sep='')

  if (length(x$n)==1)  {
    z <- sign(x$exp - x$obs) * sqrt(x$chisq)
    temp <- c(x$obs, x$exp, z, signif(1-pchisq(x$chisq, 1),digits))
    names(temp) <- c("Observed", "Expected", "Z", "p")
    print(temp)
  }
  else {
    if (is.matrix(x$obs)){
      otmp <- apply(x$obs,1,sum)
      etmp <- apply(x$exp,1,sum)
    }
    else {
      otmp <- x$obs
      etmp <- x$exp
    }
    df <- (sum(1*(etmp>0))) -1
    temp <- cbind(x$n, otmp, etmp, ((otmp-etmp)^2)/ etmp,
                  ((otmp-etmp)^2)/ diag(x$var))
    dimnames(temp) <- list(names(x$n), c("N", "Observed", "Expected",
                                         "(O-E)^2/E", "(O-E)^2/V"))
    print(temp)
    cat("\n Chisq=", format(round(x$chisq,1)),
        " on", df, "degrees of freedom, p=",
        format(signif(1-pchisq(x$chisq, df),digits)), "\n",
        "                Permutation p-value = ", x$perm.p, "\n")
  }
  invisible(x)
}

# stratified permutation function
stratperm1 <- function(n, ii, nstrat, nstrat2, uu) {
  zz <- .Fortran("strperm1",
                 as.integer(n),
                 pii=as.integer(ii),
                 as.integer(nstrat+1),
                 as.integer(nstrat2),
                 as.double(uu),
                 PACKAGE="clinfun")
  zz$pii
}


# The exact/permutation version of Jonckheere-Terpstra test
# asymptotic version is equivalent to cor.test for Kendall's tau

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

jonckheere.test <- function(x, g, alternative = c("two.sided", "less", "greater")) {
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
  STATISTIC <- jtrsum
  names(STATISTIC) <- "JT"
  if (TIES | (n > 50)) {
    STATISTIC <- (STATISTIC-jtmean)/sqrt(jtvar)
    PVAL <- 1 - pnorm(STATISTIC)
    PVAL <- switch(alternative,
                   "two.sided" = 2*min(PVAL, 1-PVAL),
                   "greater" = PVAL,
                   "less" = 1-PVAL)
  } else {
    lPVAL <- sum(djonckheere(gsize)[1:(jtrsum+1)])
    gPVAL <- 1-sum(djonckheere(gsize)[1:(jtrsum)])
    PVAL <- switch(alternative,
                   "two.sided" = 2*min(lPVAL, gPVAL),
                   "greater" = gPVAL,
                   "less" = lPVAL)
  }
  RVAL <- list(statistic = STATISTIC,
               p.value = as.numeric(PVAL),
               alternative = alternative,
               method = METHOD)
  class(RVAL) <- "htest"
  RVAL
}
