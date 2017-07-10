\name{cusum}
\alias{cusum}
\alias{cusum.qcc}
\alias{print.cusum.qcc}
\alias{summary.cusum.qcc}
\alias{plot.cusum.qcc}

\title{Cusum chart}
\description{Create an object of class 'cusum.qcc' to compute a Cusum chart for statistical quality control.}

\usage{
cusum(data, sizes, center, std.dev, head.start = 0,
      decision.interval = 5, se.shift = 1, data.name, labels,
      newdata, newsizes, newlabels, plot = TRUE, \dots)

\method{print}{cusum.qcc}(x, \dots)

\method{summary}{cusum.qcc}(object, digits = getOption("digits"), \dots)

\method{plot}{cusum.qcc}(x, add.stats = TRUE, chart.all = TRUE, 
     label.bounds = c("LDB", "UDB"), title, xlab, ylab, ylim, 
     axes.las = 0, digits =  getOption("digits"),
     restore.par = TRUE, \dots)
}
\arguments{
\item{data}{a data frame, a matrix or a vector containing observed data for the variable to chart. Each row of a data frame or a matrix, and each value of a vector, refers to a sample or ''rationale group''.}

\item{sizes}{a value or a vector of values specifying the sample sizes associated with each group. If not provided the sample sizes are obtained counting the non-\code{NA} elements of each row of a data frame or a matrix; sample sizes are set all equal to one if \code{data} is a vector.}

\item{center}{a value specifying the center of group statistics or the ''target'' value of the process.}

\item{std.dev}{a value or an available method specifying the within-group standard deviation(s) of the process. \cr
Several methods are available for estimating the standard deviation. See \code{\link{sd.xbar}} and \code{\link{sd.xbar.one}} for, respectively, the grouped data case and the individual observations case.
}

\item{head.start}{The initializing value for the above-target and
below-target cumulative sums, measured in standard errors of the summary
statistics. Use zero for the traditional Cusum chart, or a positive
value less than the \code{decision.interval} for a Fast Initial Response.}

\item{decision.interval}{A numeric value specifying the number of standard errors of the summary statistics at which the cumulative sum is out of control.}

\item{se.shift}{The amount of shift to detect in the process, measured in standard errors of the summary statistics.}

\item{data.name}{a string specifying the name of the variable which appears on the plots. If not provided is taken from the object given as data.}

\item{labels}{a character vector of labels for each group.}

\item{newdata}{a data frame, matrix or vector, as for the \code{data} argument, providing further data to plot but not included in the computations.}

\item{newsizes}{a vector as for the \code{sizes} argument providing further data sizes to plot but not included in the computations.}

\item{newlabels}{a character vector of labels for each new group defined in the argument \code{newdata}.}

\item{plot}{logical. If \code{TRUE} a Cusum chart is plotted.}

\item{add.stats}{a logical value indicating whether statistics and other information should be printed at the bottom of the chart.}

\item{chart.all}{a logical value indicating whether both statistics for \code{data} and for \code{newdata} (if given) should be plotted.}

\item{label.bounds}{a character vector specifying the labels for the the decision interval boundaries.}

\item{title}{a string giving the label for the main title.}

\item{xlab}{a string giving the label for the x-axis.}
 
\item{ylab}{a string giving the label for the y-axis.}

\item{ylim}{a numeric vector specifying the limits for the y-axis.}

\item{axes.las}{numeric in \{0,1,2,3\} specifying the style of axis labels. See \code{help(par)}.}

\item{digits}{the number of significant digits to use.}

\item{restore.par}{a logical value indicating whether the previous \code{par} settings must be restored. If you need to add points, lines, etc. to a control chart set this to \code{FALSE}.}

\item{object}{an object of class 'cusum.qcc'.}

\item{x}{an object of class 'cusum.qcc'.}

\item{\dots}{additional arguments to be passed to the generic function.} 
}
\details{Cusum charts display how the group summary statistics deviate above or below the process center or target value, relative to the standard errors of the summary statistics. Useful to detect small and permanent variation on the mean of the process.
}
\value{Returns an object of class 'cusum.qcc'.}
\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM. \cr
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons. \cr
Ryan, T. P. (2000), \emph{Statistical Methods for Quality Improvement}, 2nd ed. New York: John Wiley & Sons, Inc. \cr
Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17. \cr
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
attach(pistonrings)
diameter <- qcc.groups(diameter, sample)

q <- cusum(diameter[1:25,], decision.interval = 4, se.shift = 1)
summary(q)

q <- cusum(diameter[1:25,], newdata=diameter[26:40,])
summary(q)
plot(q, chart.all=FALSE)

detach(pistonrings)
}
\keyword{htest}
\keyword{hplot}
