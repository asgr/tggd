\name{tggd_log}
\alias{tggd_log}
\alias{dtggd_log}
\alias{ptggd_log}
\alias{qtggd_log}
\alias{rtggd_log}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Truncated Generalised Gamma Distribution
}
\description{
Density, distribution function, quantile function and random generation for the Truncated Generalised Gamma Distribution in log base 10 (log10) space.
}
\usage{
dtggd_log(x, scale=14, a=-1, b=1, xmin=10, log = FALSE)
ptggd_log(q, scale=14, a=-1, b=1, xmin=10, lower.tail = TRUE, log.p = FALSE)
qtggd_log(p, scale=14, a=-1, b=1, xmin=10, lower.tail = TRUE, log.p = FALSE, res.approx=1e-2)
rtggd_log(n, scale=14, a=-1, b=1, xmin=10, res.approx=1e-2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{
Vector of quantiles.
}
  \item{p}{
Vector of probabilities.
}
  \item{n}{
Number of observations. If length(n) > 1, the length is taken to be the number required.
}
  \item{scale}{
Vector of values for scale, which controls the transition regime between the power law slope and the exponential cut-off of the TGGD. This is analogous to the scale parameter for the standard Gamma distribution (see \code{\link{GammaDist}}).
}
  \item{a}{
Vector of values for a, which controls the power law slope of the TGGD.
}
  \item{b}{
Vector of values for b, which controls the exponential cutoff of the TGGD.
}
  \item{xmin}{
Vector of values for xmin, which controls the lower limit at which to trancate the TGGD.
}
  \item{res.approx}{
The resolution used to create the inverted CDF required to map probability integrals back onto quantiles.
}
}
\details{
This distribution function is described in detail in Murray, Robotham and Power 2016.
}
\value{
dtggd_log gives the density, ptggd_log gives the distribution function, qtggd_log gives the quantile function, and rtggd_log generates random deviates.

Invalid arguments will result in return value NaN, with a warning.

The length of the result is determined by n for rtggd_log, and is the maximum of the lengths of the numerical arguments for the other functions.

The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
}
\references{
Murray, Robotham and Power (2016)
}
\author{
Aaron Robotham
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GammaDist}} regarding the Gamma distribution.
\code{\link{RNG}} about random number generation in R.
\code{\link{Distributions}} for other standard distributions.
}
\examples{
r <- rtggd_log(20)
##Ideally the output below should equal 0, in practice it will be very close:
qtggd_log(ptggd_log(r))-r

#These should be the same:
integrate(dtggd_log,lower=10,upper=11,a=-1.5,b=0.7,xmin=10)
ptggd_log(11,a=-1.5,b=0.7,xmin=10)

#This should be very close to 1 (for a true PDF):
ptggd_log(18,a=-1.5,b=0.7,xmin=10)

\keyword{ ~tggd_log }
\keyword{ ~distribution }
