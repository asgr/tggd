\name{tggd}
\alias{tggd}
\alias{dtggd}
\alias{ptggd}
\alias{qtggd}
\alias{rtggd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Truncated Generalised Gamma Distribution
}
\description{
Density, distribution function, quantile function and random generation for the Truncated Generalised Gamma Distribution in linear space.
}
\usage{
dtggd(x, scale=1e14, a=-1, b=1, xmin=1e10, log = FALSE)
ptggd(q, scale=1e14, a=-1, b=1, xmin=1e10, lower.tail = TRUE, log.p = FALSE)
qtggd(p, scale=1e14, a=-1, b=1, xmin=1e10, lower.tail = TRUE, log.p = FALSE, res.approx=1e-2)
rtggd(n, scale=1e14, a=-1, b=1, xmin=1e10, res.approx=1e-2)
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
dtggd gives the density, ptggd gives the distribution function, qtggd gives the quantile function, and rtggd generates random deviates.

Invalid arguments will result in return value NaN, with a warning.

The length of the result is determined by n for rtggd, and is the maximum of the lengths of the numerical arguments for the other functions.

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
r <- rtggd(20)
##Ideally the output below should equal 1, in practice it will be very close:
qtggd(ptggd(r))/r

#These should be the same:
integrate(dtggd,lower=1e10,upper=1e11,a=-1.5,b=0.7,xmin=1e10)
ptggd(1e11,a=-1.5,b=0.7,xmin=1e10)

#This should be very close to 1 (for a true PDF):
ptggd(1e18,a=-1.5,b=0.7,xmin=1e10)

\keyword{ ~tggd }
\keyword{ ~distribution }
