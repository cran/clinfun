roc.curve <- function(marker, status, method=c("empirical")) {
  method <- match.arg(method)
  ux <- sort(unique(marker))
  n <- length(marker)
  n1 <- sum(status)
  n0 <- n - n1
  x0 <- marker[status == 0]
  x1 <- marker[status == 1]
  out <- NULL
  out$tpr <- c(sapply(ux, function(y, x) {
    sum(x >= y)
  } , x1)/n1, 0)
  out$fpr <- c(sapply(ux, function(y, x) {
    sum(x >= y)
  } , x0)/n0, 0)
  out$marker <- marker
  out$status <- status
  class(out) <- "roc.curve"
  out
}

print.roc.curve <- function(x, ...) {
  out <- roc.area.test(x$marker, x$status)
  cat("  ROC curve with AUC =", out$area, "and s.e. =", sqrt(out$var), "\n")
}

plot.roc.curve <- function(x, ...) {
  plot(x$fpr, x$tpr, xlab="False positive rate", ylab="True positive rate", type="l", ...)
}

lines.roc.curve <- function(x, ...) {
  lines(x$fpr, x$tpr, ...)
}
