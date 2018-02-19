\name{stats.xbar}
\alias{stats.xbar}
\alias{sd.xbar}
\alias{limits.xbar}

\title{Statistics used in computing and drawing a Shewhart xbar chart}
\description{These functions are used to compute statistics required by the xbar chart.}
\usage{
stats.xbar(data, sizes)
sd.xbar(data, sizes, 
        std.dev = c("UWAVE-R", "UWAVE-SD", "MVLUE-R", "MVLUE-SD", "RMSDF"), 
        \dots)
limits.xbar(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
}
\arguments{
  \item{data}{the observed data values}
  \item{center}{sample/group center statistic}
  \item{sizes}{samples sizes. Optional}
  \item{std.dev}{within group standard deviation. Optional for \code{sd.xbar} function, required for \code{limits.xbar}. See details.}
  \item{nsigmas}{a numeric value specifying the number of sigmas to use for computing control limits. It is ignored when the \code{conf} argument is  provided.} 
  \item{conf}{a numeric value \eqn{\in (0,1)} specifying the confidence level to use for computing control limits.}
  \item{\dots}{catches further ignored arguments.}
}
\details{
The following methods are available for estimating the process standard deviation:
\describe{
\item{\code{"UWAVE-R"}}{UnWeighted AVErage of within-group estimates based on within-group Ranges}
\item{\code{"UWAVE-SD"}}{UnWeighted AVErage of within-group estimates based on within-group Standard Deviations}
\item{\code{"MVLUE-R"}}{Minimum Variance Linear Unbiased Estimator computed as a weighted average of within-group estimates based on within-group Ranges}
\item{\code{"MVLUE-SD"}}{Minimum Variance Linear Unbiased Estimator computed as a weighted average of within-group estimates based on within-group Standard Deviations}
\item{\code{"RMSDF"}}{Root-Mean-Square estimator computed as a weighted average of within-group estimates based on within-group Standard Deviations}
}

Depending on the chart, a method may be available or not, or set as the default according to the following table:
\tabular{lcccc}{
Method \tab\tab \code{"xbar"} \tab \code{"R"} \tab \code{"S"} \cr
\code{"UWAVE-R"} \tab\tab default \tab default \tab not available \cr
\code{"UWAVE-SD"} \tab\tab  \tab not available \tab default \cr
\code{"MVLUE-R"} \tab\tab \tab \tab not available  \cr
\code{"MVLUE-SD"} \tab\tab \tab not available  \tab \cr
\code{"RMSDF"} \tab\tab \tab not available \tab
}

Detailed definitions of formulae implemented are available in the SAS/QC User's Guide.
}
\value{
The function \code{stats.xbar} returns a list with components \code{statistics} and \code{center}.

The function \code{sd.xbar} returns \code{std.dev} the standard deviation of the statistic charted. This is based on results from Burr (1969).

The function \code{limits.xbar} returns a matrix with lower and upper control limits. 
}
\references{
Burr, I.W. (1969) Control charts for measurements with varying sample sizes. \emph{Journal of Quality Technology}, 1(3), 163-167.

Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons.

Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{}
\seealso{\code{\link{qcc}} }
%\examples{}
\keyword{htest}
\keyword{hplot}

