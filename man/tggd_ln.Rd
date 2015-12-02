\name{tggd_ln}
\alias{tggd_ln}
\alias{dtggd_ln}
\alias{ptggd_ln}
\alias{qtggd_ln}
\alias{rtggd_ln}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Truncated Generalised Gamma Distribution
}
\description{
Density, distribution function, quantile function and random generation for the Truncated Generalised Gamma Distribution in natural log (ln) space. Specifically, if exp(x) is drawn from a TGGD distribution (in real space), these functions give the distribution of x, using the same parameter values.
}
\usage{
dtggd_ln(x, scale=log(1e14), a=-1, b=1, xmin=log(1e10), log = FALSE)
ptggd_ln(q, scale=log(1e14), a=-1, b=1, xmin=log(1e10), lower.tail = TRUE, log.p = FALSE)
qtggd_ln(p, scale=log(1e14), a=-1, b=1, xmin=log(1e10), lower.tail = TRUE, log.p = FALSE,
res.approx=1e-2)
rtggd_ln(n, scale=log(1e14), a=-1, b=1, xmin=log(1e10), res.approx=1e-2)
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
This distribution function is described in detail in Murray, Robotham and Power 2016.

\deqn{f(x;a,b,s,m) = \frac{b \exp((x-s)^{(a+1)}) \exp(-\exp(b(x-s)))}{s \Gamma(\frac{a+1}{b},\exp(m-s)^b)}}{f(x;a,b,s,m) = (b.exp(x-s)^(a+1).exp(-exp(b(x-s))))/gamma_inc((a+1)/b,exp(m-s)^b)}

where, from the argument list above, we use x=x, a=a, b=b, s=scale and m=xmin. \eqn{\Gamma}{gamma_inc} is the upper incomplete Gamma function as defined by the gsl \code{\link{gamma_inc}} function, using the same argument ordering, where \code{gamma_inc(a,x)==pgamma(x,a,lower=FALSE)*gamma(x)} for a>0. \code{\link{gamma_inc}} is used because it allows for the computation of upper incomplete integrals in cases where a<=0.
}
\value{
dtggd_ln gives the density, ptggd_ln gives the distribution function, qtggd_ln gives the quantile function, and rtggd_ln generates random deviates.

Invalid arguments will result in return value NaN, with a warning.

The length of the result is determined by n for rtggd_ln, and is the maximum of the lengths of the numerical arguments for the other functions.

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
r <- rtggd_ln(100,a=-2)
hist(r)

##Ideally the output below should equal 0, in practice it will be very close:
qtggd_ln(ptggd_ln(r))-r

#These should be the same:
integrate(dtggd_ln,lower=log(1e10),upper=log(1e11),a=-1.5,b=0.7,xmin=log(1e10))
ptggd_ln(log(1e11),a=-1.5,b=0.7,xmin=log(1e10))

#This should be very close to 1 (for a true PDF):
ptggd_ln(log(1e18),a=-1.5,b=0.7,xmin=log(1e10))

#To show the link to the linear and log10 (called log) variants (and the slight
#inaccuracies) these outputs should be a sequence from 0 to 1 (by=0.1):
ptggd_log(log10(qtggd(seq(0,1,by=0.1))))
ptggd_log(qtggd_ln(seq(0,1,by=0.1))/log(10))
}
\keyword{ ~tggd_ln }
\keyword{ ~distribution }
