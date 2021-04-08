\name{ewma}
\alias{ewma}
\alias{ewma.qcc}
\alias{print.ewma.qcc}
\alias{summary.ewma.qcc}
\alias{plot.ewma.qcc}

\title{EWMA chart}
\description{Create an object of class \code{'ewma.qcc'} to compute and draw an Exponential Weighted Moving Average (EWMA) chart for statistical quality control.}

\usage{
ewma(data, sizes, center, std.dev, 
     lambda = 0.2, nsigmas = 3, 
     newdata, newsizes, \dots)

\method{print}{ewma.qcc}(x, digits = getOption("digits"), \dots)

\method{plot}{ewma.qcc}(x, 
     xtime, 
     add.stats = qcc.options("add.stats"), 
     chart.all = qcc.options("chart.all"), 
     fill = qcc.options("fill"),
     label.center = "CL",
     label.limits = c("LCL", "UCL"), 
     title, xlab, ylab, xlim, ylim, 
     digits = getOption("digits"), \dots)
}

\arguments{
\item{data}{a data frame, a matrix or a vector containing observed data for the variable to chart. Each row of a data frame or a matrix, and each value of a vector, refers to a sample or ''rationale group''.}

\item{sizes}{a value or a vector of values specifying the sample sizes associated with each group. If not provided the sample sizes are obtained counting the non-\code{NA} elements of each row of a data frame or a matrix; sample sizes are set all equal to one if \code{data} is a vector.}

\item{center}{a value specifying the center of group statistics or target.}

\item{std.dev}{a value or an available method specifying the within-group standard deviation(s) of the process. \cr
Several methods are available for estimating the standard deviation. See \code{\link{sd.xbar}} and \code{\link{sd.xbar.one}} for, respectively, the grouped data case and the individual observations case.
}

\item{lambda}{the smoothing parameter \eqn{0 \le \lambda \le 1}{0 <= lambda <= 1}}

\item{nsigmas}{a numeric value specifying the number of sigmas to use for computing control limits.}

\item{newdata}{a data frame, matrix or vector, as for the \code{data} argument, providing further data to plot but not included in the computations.}

\item{newsizes}{a vector as for the \code{sizes} argument providing further data sizes to plot but not included in the computations.}

\item{xtime}{a vector of date-time values as returned by \code{\link{Sys.time}} and \code{\link{Sys.Date}}. If provided it is used for x-axis so it must be of the same length as the statistic charted.}

\item{add.stats}{a logical value indicating whether statistics and other information should be printed at the bottom of the chart.}

\item{chart.all}{a logical value indicating whether both statistics for \code{data} and for \code{newdata} (if given) should be plotted.}

\item{fill}{a logical value specifying if the in-control area should be filled with the color specified in 
\code{qcc.options("zones")$fill}.}

\item{label.center}{a character  specifying the label for center line.}

\item{label.limits}{a character vector specifying the labels for control limits.}

\item{title}{a character string specifying the main title. Set \code{title = NULL} to remove the title.}

\item{xlab, ylab}{a string giving the label for the x-axis and the y-axis.}

\item{xlim, ylim}{a numeric vector specifying the limits for the x-axis and the y-axis.}

\item{digits}{the number of significant digits to use.}

\item{x}{an object of class {'ewma.qcc'}.}

\item{\dots}{additional arguments to be passed to the generic function.} 
}

\details{EWMA chart smooths a series of data based on a moving average with weights which decay exponentially. Useful to detect small and permanent variation on the mean of the process.
}
\value{Returns an object of class {'ewma.qcc'}.}

\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM.

Montgomery, D.C. (2013) \emph{Introduction to Statistical Quality Control}, 7th ed. New York: John Wiley & Sons.

Ryan, T. P. (2011), \emph{Statistical Methods for Quality Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.

Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17.

Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{ ~~further notes~~ }
\seealso{\code{\link{qcc}}, \code{\link{ewmaSmooth}}, \code{\link{cusum}}}
\examples{
##
## Grouped-data
##
data(pistonrings)
diameter <- qccGroups(data = pistonrings, diameter, sample)

q = ewma(diameter[1:25,], lambda=0.2, nsigmas=3)
summary(q)
plot(q)

ewma(diameter[1:25,], lambda=0.2, nsigmas=2.7, newdata=diameter[26:40,]) 

##
## Individual observations
##
data(viscosity)
with(viscosity, 
     ewma(viscosity[trial], lambda = 0.2, nsigmas = 2.7,
          newdata = viscosity[!trial]))
}
\keyword{htest}
\keyword{hplot}
