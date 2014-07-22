\name{gsdesign-internal}
\title{Group Sequential Designs Internal Functions}
\alias{gsd.drift}
\description{
 Internal gsdesign functions.
}
\usage{
gsd.bdryconstant(ifrac, eprob = 0.05, delta = 0.5,
    alternative = c("two.sided", "one.sided"), tol = 0.00001, ...)
gsd.drift.efficacy(ifrac, delta.eb, sig.level=0.05, pow = 0.8,
    alternative = c("two.sided", "one.sided"), tol = 0.00001, ...)
gsd.drift.both(ifrac, delta.eb, delta.fb, sig.level = 0.05, pow=0.8,
    alternative = c("two.sided", "one.sided"), tol = 0.00001, ...)
gsd.drift(ifrac, sig.level = 0.05, pow = 0.8, delta.eb = 0.5, delta.fb = NULL,
    alternative = c("two.sided", "one.sided"), tol = 0.00001, ...)
}
\details{
  These functions are not to be called by the user.
}
\keyword{internal}
