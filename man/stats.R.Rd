\name{stats.R}
\alias{stats.R}
\alias{sd.R}
\alias{limits.R}

\title{Statistics used in computing and drawing a Shewhart R chart}
\description{These functions are used to compute statistics required by the R chart.}
\usage{
stats.R(data, sizes)
sd.R(data, sizes, std.dev = c("UWAVE-R", "MVLUE-R"), \dots)
limits.R(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
}
\arguments{
  \item{data}{the observed data values}
  \item{center}{sample/group center statistic.}
  \item{sizes}{samples sizes. Optional}
  \item{std.dev}{within group standard deviation. Optional for \code{sd.R} function, required for \code{limits.R}. See \code{\link{sd.xbar}}.}
  \item{nsigmas}{a numeric value specifying the number of sigmas to use for computing control limits. It is ignored when the \code{conf} argument is  provided.} 
  \item{conf}{a numeric value \eqn{\in (0,1)} specifying the confidence level to use for computing control limits.}
  \item{\dots}{catches further ignored arguments.}
}
%\details{}
\value{
The function \code{stats.R} returns a list with components \code{statistics} and \code{center}.

The function \code{sd.R} returns \code{std.dev} the standard deviation of the statistic charted.

The function \code{limits.R} returns a matrix with lower and upper control limits. 
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

