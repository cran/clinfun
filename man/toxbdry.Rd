\name{toxbdry}
\title{Stopping rule for toxicity monitoring}
\alias{toxbdry}
\alias{bdrycross.prob}
\alias{print.toxbdry}
\keyword{design}
\description{
  Computes a stopping rule and its operating characteristics for
  toxicity monitoring based repeated significance testing.
}
\usage{
toxbdry(pLo, pHi, n, cP0=0.1, cP1=0.9, ngrid=6, niter=5)
bdrycross.prob(n, r, ptox)
\method{print}{toxbdry}(x, \dots)
}
\arguments{
  \item{pLo}{the toxicity rate that is acceptable.}
  \item{pHi}{the toxicity rate that is too high and hence unacceptable.}
  \item{n}{vector of times (sample size) when toxicty is moniroted.}
  \item{r}{vector of maximum acceptable toxicities corresponding to n.}
  \item{ptox}{the toxicity rates for which the operating characteristics
    are calculated.}
  \item{cP0}{boundary crossing probability under pLo i.e. type I error
    or the probability of declaring a treatment with toxicity rate pLo
    unacceptable.}
  \item{cP1}{boundary crossing probability under pHi i.e. power or the
    probability of declaring a treatment with toxicity rate pHi
    unacceptable.}
  \item{ngrid}{the number of toxicity rates from pLo to pHi for which
    the operating characteristics are computed.}
  \item{niter}{the number of iterations run to obtain the boundary.}
  \item{x}{object returned by the function toxbdry.}
  \item{...}{additional arguments to print.}
}
\value{
  the function returns a list with:
  \item{looks}{when toxicty is monitored - same as input n.}
  \item{lo.bdry}{lower boundary is a vector of maximum acceptable number
    of toxicities corresponding the number of subjects in n.  The
    boundary crossing probability for this is slightly above cP0.}
  \item{hi.bdry}{upper boundary is a vector of maximum acceptable number
    of toxicities corresponding the number of subjects in n.  The
    boundary crossing probability for this is slightly below cP0.}
  \item{bdry.oc}{the operating characteristics i.e the toxicity
    rate, the probability of crossing, stopping (i.e. cross before the
    last observation) and the expected sample size for both the low (lo)
    and high (hi) boundaries.}
  \item{bdry.alpha}{the alpha levels for testing at each look for the
    two boundaries.}

  stopping for toxicity is done when the number of toxicities exceeda
  the boundary i.e. the boundary gives the maximum acceptable number.
}
\examples{
  toxbdry(0.2, 0.35, c(20,40,60,75))
  toxbdry(0.2, 0.3, c(20,40,60,75), cP0=0.15, cP1=0.8)
  toxbdry(0.1, 0.3, 1:30) # continuous monitoring
}
