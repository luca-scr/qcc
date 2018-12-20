\name{stats.xbar.one}
\alias{stats.xbar.one}
\alias{sd.xbar.one}
\alias{limits.xbar.one}

\title{Statistics used in computing and drawing a Shewhart xbar chart for one-at-time data}
\description{These functions are used to compute statistics required by the xbar chart for one-at-time data.}
\usage{
stats.xbar.one(data, sizes)
sd.xbar.one(data, sizes, std.dev = c("MR", "SD"), r = 2, \dots)
limits.xbar.one(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
}
\arguments{
  \item{data}{the observed data values}
  \item{center}{sample/group center statistic.}
  \item{sizes}{samples sizes. Not needed, \code{size = 1} is used.}
  \item{r}{number of successive pairs of observations for computing the standard deviation based on moving ranges of r points.}
  \item{std.dev}{within group standard deviation. Optional for \code{sd.xbar.one} function, required for \code{limits.xbar.one}. See details.}
   \item{nsigmas}{a numeric value specifying the number of sigmas to use for computing control limits. It is ignored when the \code{conf} argument is  provided.} 
  \item{conf}{a numeric value in \eqn{(0,1)} specifying the confidence level to use for computing control limits.}
  \item{\dots}{catches further ignored arguments.}
}
\details{
Methods  available for estimating the process standard deviation:
\itemize{
\item{\code{"MR"} = moving range:} this is estimate is based on the scaled mean of moving ranges
\item{\code{"SD"} = sample standard deviation:} this estimate is defined as \code{sd(x)/cd(n)}, where \code{n} is the  number of individual measurements of \code{x}.
}
}
\value{
The function \code{stats.xbar.one} returns a list with components \code{statistics} and \code{center}.

The function \code{sd.xbar.one} returns \code{std.dev} the standard deviation of the statistic charted.

The function \code{limits.xbar.one} returns a matrix with lower and upper control limits. 
}
\references{
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons.

Ryan, T. P. (2011), \emph{Statistical Methods for Quality Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.

Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{}
\seealso{\code{\link{qcc}} }
\examples{
# Water content of antifreeze data (Wetherill and Brown, 1991, p. 120)
x <- c(2.23, 2.53, 2.62, 2.63, 2.58, 2.44, 2.49, 2.34, 2.95, 2.54, 2.60, 2.45,
       2.17, 2.58, 2.57, 2.44, 2.38, 2.23, 2.23, 2.54, 2.66, 2.84, 2.81, 2.39,
       2.56, 2.70, 3.00, 2.81, 2.77, 2.89, 2.54, 2.98, 2.35, 2.53)
# the Shewhart control chart for one-at-time data
# 1) using MR (default)
qcc(x, type="xbar.one", data.name="Water content (in ppm) of batches of antifreeze")
# 2) using SD
qcc(x, type="xbar.one", std.dev = "SD", data.name="Water content (in ppm) of batches of antifreeze")

# "as the size increases further, we would expect sigma-hat to settle down
#  at a value close to the overall sigma-hat" (Wetherill and Brown, 1991,
# p. 121)
sigma <- NA
k <- 2:24
for (j in k)
    sigma[j] <- sd.xbar.one(x, k=j)
plot(k, sigma[k], type="b")     # plot estimates of sigma for 
abline(h=sd(x), col=2, lty=2)   # different values of k
}
\keyword{htest}
\keyword{hplot}

