\name{stats.c}
\alias{stats.c}
\alias{sd.c}
\alias{limits.c}

\title{Functions to plot Shewhart c chart}
\description{Statistics used in computing and drawing a Shewhart c chart.}
\usage{
stats.c(data, sizes)
sd.c(data, sizes, ...)
limits.c(center, std.dev, sizes, conf)
}
\arguments{
  \item{data}{the observed data values}
  \item{center}{sample/group center statistic.}
  \item{sizes}{samples sizes.}
  \item{std.dev}{within group standard deviation.}
  \item{conf}{a numeric value used to compute control limits, specifying the number of standard deviations (if \code{conf} > 1) or the confidence level (if 0 < \code{conf} < 1).}
  \item{...}{further arguments are ignored.}
}
%\details{}
\value{
The function \code{stats.c} returns a list with components \code{statistics} and \code{center}.

The function \code{sd.c} returns \code{std.dev} the standard deviation of the statistic charted.

The function \code{limits.c} returns a matrix with lower and upper control limits. 
}
\references{
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{}
\seealso{\code{\link{qcc}} }
%\examples{}
\keyword{htest}
\keyword{hplot}

