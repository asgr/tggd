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
Density, distribution function, quantile function and random generation for the Truncated Generalised Gamma Distribution in log base 10 (log10) space. Specifically, if 10^x is drawn from a TGGD distribution (in real space), these functions give the distribution of x, using the same parameter values.
}
\usage{
dtggd_log(x, scale=14, a=-1, b=1, xmin=10, log = FALSE)
ptggd_log(q, scale=14, a=-1, b=1, xmin=10, lower.tail = TRUE, log.p = FALSE)
qtggd_log(p, scale=14, a=-1, b=1, xmin=10, lower.tail = TRUE, log.p = FALSE,
res.approx=1e-2)
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
  \item{log, log.p}{
Logical; if TRUE, probabilities/densities p are returned as log(p).
}
  \item{lower.tail}{
Logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
}
}
\details{
This distribution function is described in detail in Murray, Robotham and Power 2016. The PDF is given by:

\deqn{f(x;a,b,s,m) = \frac{\ln(10).b (10^{(x-s)})^{(a+1)} \exp(-10^{(b(x-s))})}{s \Gamma(\frac{a+1}{b},(10^{(m-s)})^b)}}{f(x;a,b,s,m) = (ln(10).b.(10^(x-scale))^(a+1).exp(-10^(b(x-s))))/gamma_inc((a+1)/b,(10^(m-s))^b)}

where, from the argument list above, we use x=x, a=a, b=b, s=scale and m=xmin. \eqn{\Gamma}{gamma_inc} is the upper incomplete Gamma function as defined by the gsl \code{\link{gamma_inc}} function, using the same argument ordering, where \code{gamma_inc(a,x)==pgamma(x,a,lower=FALSE)*gamma(x)} for a>0. \code{\link{gamma_inc}} is used because it allows for the computation of upper incomplete integrals in cases where a<=0.
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
Aaron Robotham, Steven Murray
}
\note{
The intended application of the log-space version of the TGGD is to provide the correct distribution when variates are drawn from a real-space TGGD, but there are priors on their uncertainty which operate in logspace (eg. a lognormal distribution). The likelihood of a given set of parameters is incorrect in such a case if the real-space version is used without an adjustment to the Jacobian.
}
\seealso{
\code{\link{GammaDist}} regarding the Gamma distribution.
\code{\link{RNG}} about random number generation in R.
\code{\link{Distributions}} for other standard distributions.
}
\examples{
r <- rtggd_log(100,a=-2)
hist(r)

##Ideally the output below should equal 0, in practice it will be very close:
qtggd_log(ptggd_log(r))-r

#These should be the same:
integrate(dtggd_log,lower=10,upper=11,a=-1.5,b=0.7,xmin=10)
ptggd_log(11,a=-1.5,b=0.7,xmin=10)

#This should be very close to 1 (for a true PDF):
ptggd_log(18,a=-1.5,b=0.7,xmin=10)

#To show the link to the linear and ln variants (and the slight inaccuracies) these
#outputs should be a sequence from 0 to 1 (by=0.1):
ptggd_ln(log(qtggd(seq(0,1,by=0.1))))
ptggd_ln(qtggd_log(seq(0,1,by=0.1))*log(10))
}
\keyword{ ~tggd_log }
\keyword{ ~distribution }
