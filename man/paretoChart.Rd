\name{paretoChart}
\alias{paretoChart}
\alias{print.paretoChart}
\alias{plot.paretoChart}
\alias{blues.colors}

\title{Pareto chart}
\description{Computes a table of statistics and plot a Pareto chart.}

\usage{
paretoChart(data, \dots)

\method{plot}{paretoChart}(x, title, xlab, 
     ylab = "Frequency", ylab2 = "Cumulative percentage", 
     ylim, col = blues.colors(nlevels), \dots) 
}
\arguments{
\item{data}{a vector of values. \code{names(data)} are used for labelling the bars.}
\item{plot}{a logical specifying if the chart should be provided (\code{TRUE}, default).}
\item{x}{an object of class \code{'paretoChart'} returned by a call to \code{paretoChart()} function.}
\item{title}{a character string specifying the main title. Set \code{title = NULL} to remove the title.}
\item{xlab}{a string specifying the label for the x-axis.}
\item{ylab}{a string specifying the label for the y-axis.}
\item{ylab2}{a string specifying the label for the second y-axis on the right side.}
\item{ylim}{a numeric vector specifying the limits for the y-axis.}
\item{col}{a value for the color, a vector of colors, or a palette for the bars. See the help for \code{\link{colors}} and \code{\link{palette}}.}
\item{\dots}{catch other optional arguments.}
}

\value{Returns an object of class \code{'paretoChart'} containing the descriptive statistics used to draw the Pareto chart. This object has associated a \code{print} and \code{plot} method.}

\details{A Pareto chart is a barplot where the categories are ordered in non increasing order, and a line is also added to show the cumulative sum.}

\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM.

Montgomery, D.C. (2013) \emph{Introduction to Statistical Quality Control}, 7th ed. New York: John Wiley & Sons.

Ryan, T. P. (2011), \emph{Statistical Methods for Quality Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.

Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17.

Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
\seealso{\code{\link{barplot}}}
\examples{
defect <- c(80, 27, 66, 94, 33)
names(defect) <- c("price code", "schedule date", "supplier code", "contact num.", "part num.")
pc = paretoChart(defect, ylab = "Error frequency")
pc
plot(pc)

plot(paretoChart(defect, ylab = "Error frequency"), col=rainbow(length(defect)))
}
\keyword{hplot}
