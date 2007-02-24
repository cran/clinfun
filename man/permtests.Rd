\name{permtests}
\title{Permutation versions of some common tests}
\alias{permtests}
\alias{permlogrank}
\description{
  Small sample tests using the permutation reference distributions.
}
\usage{
  permlogrank(formula, data, subset, na.action, rho=0, nperm=5000)
}
\arguments{
  \item{nperm}{number of permutations for the reference distribution}
  \item{formula, data, subset, na.action, rho}{see survdiff for details}
}
\details{
  permlogrank is the permutation version of k-sample survdiff
}
\keyword{htest}
