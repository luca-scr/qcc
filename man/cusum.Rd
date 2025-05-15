\name{cusum}
\alias{cusum}
\alias{cusum.qcc}
\alias{print.cusum.qcc}
\alias{summary.cusum.qcc}
\alias{plot.cusum.qcc}

\title{Cusum chart}
\description{Create an object of class \code{'cusum.qcc'} to compute a Cusum chart for statistical quality control.}

\usage{
cusum(data, sizes, center, std.dev, 
      decision.interval = 5, se.shift = 1,
      head.start = 0,
      newdata, newsizes, \dots)

\method{print}{cusum.qcc}(x, digits = getOption("digits"), \dots)

\method{plot}{cusum.qcc}(x, 
     xtime = NULL,
     add.stats = qcc.options("add.stats"), 
     chart.all = qcc.options("chart.all"), 
     fill = qcc.options("fill"),
     label.bounds = c("LDB", "UDB"),
     title, xlab, ylab, xlim, ylim, 
     digits =  getOption("digits"), \dots)
}
\arguments{
\item{data}{a data frame, a matrix or a vector containing observed data for the variable to chart. Each row of a data frame or a matrix, and each value of a vector, refers to a sample or ''rationale group''.}

\item{sizes}{a value or a vector of values specifying the sample sizes associated with each group. If not provided the sample sizes are obtained counting the non-\code{NA} elements of each row of a data frame or a matrix; sample sizes are set all equal to one if \code{data} is a vector.}

\item{center}{a value specifying the center of group statistics or the ''target'' value of the process.}

\item{std.dev}{a value or an available method specifying the within-group standard deviation(s) of the process. \cr
Several methods are available for estimating the standard deviation. See \code{\link{sd.xbar}} and \code{\link{sd.xbar.one}} for, respectively, the grouped data case and the individual observations case.
}

\item{decision.interval}{A numeric value specifying the number of standard errors of the summary statistics at which the cumulative sum is out of control.}

\item{se.shift}{The amount of shift to detect in the process, measured in standard errors of the summary statistics.}

\item{head.start}{The initializing value for the above-target and
below-target cumulative sums, measured in standard errors of the summary
statistics. Use zero for the traditional Cusum chart, or a positive
value less than the \code{decision.interval} for a Fast Initial Response.}

\item{newdata}{a data frame, matrix or vector, as for the \code{data} argument, providing further data to plot but not included in the computations.}

\item{newsizes}{a vector as for the \code{sizes} argument providing further data sizes to plot but not included in the computations.}

\item{xtime}{a vector of date-time values as returned by \code{\link{Sys.time}} and \code{\link{Sys.Date}}. If provided it is used for x-axis so it must be of the same length as the statistic charted.}

\item{add.stats}{a logical value indicating whether statistics and other information should be printed at the bottom of the chart.}

\item{chart.all}{a logical value indicating whether both statistics for \code{data} and for \code{newdata} (if given) should be plotted.}

\item{fill}{a logical value specifying if the in-control area should be filled with the color specified in \code{qcc.options("zones")$fill}.}

\item{label.bounds}{a character vector specifying the labels for the the decision interval boundaries.}

\item{title}{a character string specifying the main title. Set \code{title = NULL} to remove the title.}

\item{xlab, ylab}{a string giving the label for the x-axis and the y-axis.}

\item{xlim, ylim}{a numeric vector specifying the limits for the x-axis and the y-axis.}

\item{digits}{the number of significant digits to use.}

\item{x}{an object of class \code{'cusum.qcc'}.}

\item{\dots}{additional arguments to be passed to the generic function.} 
}
\details{Cusum charts display how the group summary statistics deviate above or below the process center or target value, relative to the standard errors of the summary statistics. Useful to detect small and permanent variation on the mean of the process.
}
\value{Returns an object of class \code{'cusum.qcc'}.}
\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM.

Montgomery, D.C. (2013) \emph{Introduction to Statistical Quality Control}, 7th ed. New York: John Wiley & Sons.

Ryan, T. P. (2011), \emph{Statistical Methods for Quality Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.

Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17.

Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{ ~~further notes~~ }
\seealso{\code{\link{qcc}}, \code{\link{ewma}}}
\examples{
##
## Grouped-data
##
data(pistonrings)
diameter  = qccGroups(data = pistonrings, diameter, sample)

q  = cusum(diameter[1:25,], decision.interval = 4, se.shift = 1)
summary(q)
plot(q)

q  = cusum(diameter[1:25,], newdata=diameter[26:40,])
summary(q)
plot(q, chart.all=FALSE)

##
## Individual observations
##
data(viscosity)
q  = with(viscosity, cusum(viscosity[trial], newdata = viscosity[!trial]))
summary(q)
plot(q)
}
\keyword{htest}
\keyword{hplot}


