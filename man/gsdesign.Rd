\name{gsdesign}
\title{Group Sequential Designs}
\alias{gsdesign}
\alias{gsdesign.binomial}
\alias{gsdesign.normal}
\alias{gsdesign.survival}
\description{
  Functions to calculate sample size for group sequential designs
}
\usage{
gsdesign.binomial(ifrac, pC, pE, sig.level = 0.05, power = 0.8,
  delta.eb=0.5, delta.fb = 0, alternative = c("two.sided", "one.sided"),
  tol=0.001) 
gsdesign.normal(ifrac, delta, sd = 1, sig.level = 0.05, power = 0.8,
  delta.eb = 0.5, delta.fb = 0, alternative = c("two.sided", "one.sided"),
  tol=0.001)
gsdesign.survival(ifrac, haz.ratio, sig.level = 0.05, power = 0.8,
  delta.eb = 0.5, delta.fb = 0, alternative = c("two.sided", "one.sided"),
  tol=0.001)
}
\arguments{
  \item{ifrac}{information fraction ratio of current to total sample
    size of number of events.}
  \item{pC}{prob of success of the standard therapy (for binomial data)}
  \item{pE}{prob of success of the experimental therapy (for binomial data)}
  \item{delta}{true difference in means (for normal data)}
  \item{sd}{standard deviation (for normal data)}
  \item{haz.ratio}{hazard ratio (for survival comarison)}
  \item{sig.level}{significance level (type I error probability)}
  \item{power}{power of test (1 minus type II error probability)}
  \item{delta.eb}{power for efficacy boundary in the Pocock (power=0) to
    O'Brien-Fleming (power=0.5) family (default is 0.5)}
  \item{delta.fb}{power for futility boundary in the Pocock (power=0) to
    O'Brien-Fleming (power=0.5) family (default is 0.5)}
  \item{alternative}{one- or two-sided test.}
  \item{tol}{tolerance level for multivariate normal probability
    computation.}
}
\value{
  a list with ifrac, sig.level, power, alternative, delta.eb and:
  \item{efbdry}{the critical value to use at the different looks.}
  \item{sample.size}{the sample size per arm for binomial/normal data.}
  \item{num.events}{the total number of failures which should be
    converted to number of subjects using censoring proportion.}  
}  
\keyword{design}
