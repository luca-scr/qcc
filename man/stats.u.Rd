\name{stats.u}
\alias{stats.u}
\alias{sd.u}
\alias{limits.u}

\title{Statistics used in computing and drawing a Shewhart u chart}
\description{These functions are used to compute statistics required by the u chart.}
\usage{
stats.u(data, sizes)
sd.u(data, sizes, \dots)
limits.u(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
}
\arguments{
  \item{data}{the observed data values}
  \item{center}{sample/group center statistic.}
  \item{sizes}{samples sizes.}
  \item{std.dev}{within group standard deviation.}
  \item{nsigmas}{a numeric value specifying the number of sigmas to use for computing control limits. It is ignored when the \code{conf} argument is  provided.} 
  \item{conf}{a numeric value in \eqn{(0,1)} specifying the confidence level to use for computing control limits.}
  \item{\dots}{catches further ignored arguments.}
}
%\details{}
\value{
The function \code{stats.u} returns a list with components \code{statistics} and \code{center}.

The function \code{sd.u} returns \code{std.dev} the standard deviation of the statistic charted.

The function \code{limits.u} returns a matrix with lower and upper control limits. 
}
\references{
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons.

Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{}
\seealso{\code{\link{qcc}} }
%\examples{}
\keyword{htest}
\keyword{hplot}

