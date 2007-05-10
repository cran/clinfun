\name{permtests}
\title{Permutation versions of some common tests}
\alias{permtests}
\alias{permlogrank}
\alias{jonckheere.test}
\description{
  Small sample tests using the permutation reference distributions.
}
\usage{
  permlogrank(formula, data, subset, na.action, rho=0, nperm=5000)
  jonckheere.test(x, g, alternative = c("two.sided", "increasing",
                 "decreasing"))
}
\arguments{
  \item{nperm}{number of permutations for the reference distribution}
  \item{formula, data, subset, na.action, rho}{see survdiff for details}
  \item{x, g}{data and group vector}
  \item{alternative}{means are monotonic (two.sided), increasing, or
    decreasing}
}
\details{
  permlogrank is the permutation version of k-sample survdiff
  
  jonckheere.test is the exact (permutation) version of the
  Jonckheere-Terpstra test.  The asymptotic version is equivalent to
  cor.test(x, g, method="k").
}
\keyword{htest}
