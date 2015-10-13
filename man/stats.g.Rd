\name{stats.g}
\alias{stats.g}
\alias{sd.g}
\alias{limits.g}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Statistics used in computing and drawing a Shewhart g chart}
\description{
  These functions are used to compute statistics required by the g chart
  (geometric distribution) for use with the qcc package.
}
\usage{
stats.g(data, sizes)
sd.g(data, sizes, ...)
limits.g(center, std.dev, sizes, conf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ the observed data values }
  \item{center}{ sample center statistic }
  \item{sizes}{ sample sizes (not used) }
  \item{std.dev}{ standard deviation of geometric distribution }
  \item{conf}{ a numeric value used to compute control limits,
    specifying the number of standard deviations (if 'conf' > 1) or the
    confidence level (if 0 < 'conf' < 1).}
  \item{...}{further arguments are ignored.}
}
\value{
  The function \code{stats.g()} returns a list with components \code{statistics} and \code{center}.

  The function \code{sd.g()} returns \code{std.dev} the standard deviation
  \eqn{sqrt(1-p)/p}.

  The function \code{limits.g()} returns a matrix with lower and upper control limits.
}
\details{
  The g chart plots the number of non-events between events.  np charts
  do not work well when the probability of an event is rare (see example
  below).  Instead of plotting the number of events, the g chart plots
  the number of non-events between events.
}

\references{
Kaminsky, FC et. al. (1992) \emph{Statistical Control Charts Based on a Geometric Distribution}, Journal of Quality Technology, 24, pp 63--69. \cr
Yang, Z et. al. (2002) On the Performance of Geometric Charts with
Estimated Control Limits, \emph{Journal of Quality Technology}, 34, pp 448--458.
}
\author{ Greg Snow \email{greg.snow@ihc.com} }
\note{ The geometric distribution is quite skewed so it is best to set
  conf at the required confidence interval (0 < conf < 1) rather than as
a multiplier of sigma.}

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{qcc}}
\examples{
success <- rbinom(1000, 1, 0.01)
num.noevent <- diff(which(c(1,success)==1))-1
qcc(success, type = "np", sizes = 1)
qcc(num.noevent, type = "g")
}
\keyword{ hplot }% __ONLY ONE__ keyword per line
