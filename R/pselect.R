pselect <- function(n, p, min.diff=1, min.resp=0) {
  ntrt <- length(p)
  if (ntrt <= 1) stop("there should be at least 2 treatments")
  if (min(p) < 0 | max(p) > 1) stop("p should be a vector of probabilities")
  if ((min.diff != round(min.diff)) | (min.diff < 1)) stop("min.diff should be a positive integer")
  psel0 <- prod(pbinom(min.resp-1, n, p))
  psel <- rep(0, ntrt)
  for(i in max(min.diff,min.resp):n) {
    pb <- pbinom(i-min.diff, n, p)
    psel <- psel + dbinom(i,n,p)*prod(pb)/pb
  }
  out <- list()
  if (min.resp > 0)  out$prob.none.selected <- psel0
  out$prob.not.unique <- 1-sum(psel)-psel0
  out$prob.selection <- cbind(p,psel)
  colnames(out$prob.selection) <- c("resp.rate", "prob.selection")
  rownames(out$prob.selection) <- paste("trt", 1:ntrt)
  out
}
