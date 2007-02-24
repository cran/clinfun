\name{pselect}
\title{Probability of selection under pick the winner rule}
\alias{pselect}
\keyword{design}
\description{
  Calculates the probability of selecting the treatment with the higher
  response rate under the pick the winner rule.
}
\usage{
pselect(n, p, min.diff=1, min.resp=0)
}
\arguments{
  \item{n}{sample size for each treatment arm.} 
  \item{p}{vector of response rates for the treatments.}
  \item{min.diff}{the minimum number by which the number of responses
    should exceed for the treatment arm to be chosen.  This must be a
    positive integer (default=1).}
  \item{min.resp}{the minimum number of responses in each treatment arm
    for it to be considered further (default=0).}
}
\value{
  the function returns a list with:
  \item{prob.none.selected}{is the probability that no treatment has at
    least min.resp responses. this element is present only if min.resp
    is greater than 0.}
  \item{prob.not.unique}{this is the probability that the best treatment
    has min.resp responses but exceeds the second best by less than
    min.diff responses.}
  \item{prob.selection}{this is a matrix which for each treatment gives
    the response probability and the probability of selecting it
    i.e. the number of responses in the chosen arm is at least min.resp
    and exceeds that in the remaining arms by at least min.diff.}
}
\examples{
  pselect(18, c(0.2, 0.2, 0.2)) # selection when no diffrence i.e. type I error
  pselect(18, c(0.2, 0.2, 0.4)) # selection probability
  pselect(26, c(0.2, 0.2, 0.4), min.diff=2, min.resp=3)
}
\references{
  Simon R, Wittes RE, Ellenberg SS. (1985). Randomized phase II clinical
  trials. \emph{Cancer Treat Rep} 69, 1375-1381.
}
