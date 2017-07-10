\name{qcc.options}
\alias{qcc.options}

\title{Set or return options for the 'qcc' package.}

\description{This function can be used to control the behavior of the
'qcc' library such as the background color, out-of-control points appearance,
 and many others.}
\usage{qcc.options(\dots)}
\arguments{
\item{\dots}{the option to be set or retrieved. See details.}
}
\details{
The available options are:
\describe{

\item{\code{exp.R.unscaled}}{a vector specifying, for each sample size, the expected value of the relative range (i.e. \eqn{R/\sigma}) for a normal distribution. This appears as \eqn{d_2} on most tables containing factors for the construction of control charts.}

\item{\code{se.R.unscaled}}{a vector specifying, for each sample size, the standard error of the relative range (i.e. \eqn{R/\sigma}) for a normal distribution. This appears as \eqn{d_3} on most tables containing factors for the construction of control charts.}

\item{\code{beyond.limits$pch}}{plotting character used to highlight points beyond control limits.}  

\item{\code{beyond.limits$col}}{color used to highlight points beyond control limits.}  

\item{\code{violating.runs$pch}}{plotting character used to highlight points violating runs.}

\item{\code{violating.runs$col}}{color used to highlight points violating runs.}  

\item{\code{run.length}}{the maximum value of a run before to signal a point as out of control.}  

\item{\code{bg.margin}}{background color used to draw the margin of the charts.} 

\item{\code{bg.figure}}{background color used to draw the figure of the charts.}

\item{\code{cex}}{character expansion used to draw plot annotations (labels, title, tickmarks, etc.).}

\item{\code{font.stats}}{font used to draw text at the bottom of control charts.}

\item{\code{cex.stats}}{character expansion used to draw text at the bottom of control charts.}

}
}
\value{
If the functions is called with no argument return a list of available options. 

If an option argument is provided the corresponding value is returned. 

If a value is associated with an option argument, such option is set and the list of updated option values is invisibly returned.
In this case the list \code{.qcc.options} is modified and any modification will remain in effect for the rest of the session.
}

\author{Luca Scrucca}
%\note{}
\seealso{\code{\link{qcc}} }
\examples{
old <- qcc.options()			# save defaults
qcc.options("cex.stats")		# get a single parameter
qcc.options("cex.stats"=1.2)	# change parameters
qcc.options(bg.margin="azure2")
qcc.options("violating.runs" = list(pch = 15, col = "purple"))
qcc.options("beyond.limits" = list(pch = 15, col = "orangered"))
qcc(rnorm(100), type = "xbar.one", std.dev = 0.7)	# see the results
qcc.options(old)				# restore old defaults 
}
\keyword{htest}
\keyword{hplot}

