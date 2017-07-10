\name{ellipseChart}
\alias{ellipseChart}

\title{Multivariate Quality Control Charts}
\description{Plot an ellipse chart for a bivariate quality control data.}

\usage{
ellipseChart(object, chart.all = TRUE, show.id = FALSE, ngrid = 50,
             confidence.level, correct.multiple = TRUE, 
             title, xlim, ylim, xlab, ylab,
             restore.par = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class 'mqcc'.}
  \item{chart.all}{a logical value indicating whether both statistics for \code{data} and for \code{newdata} (if given) should be plotted.}
  \item{show.id}{a logical value indicating whether to plot point labels (\code{TRUE}) or symbols (\code{FALSE}) for group means.}
  \item{ngrid}{a value for the size of the grid over which the ellipse is evaluated.}
  \item{confidence.level}{a numeric value between 0 and 1 specifying the confidence level of the computed probability limits.}
  \item{correct.multiple}{a logical value indicating whether to correct or not for multiple comparisons.}
  \item{title}{a string giving the label for the main title.}
  \item{xlim}{a numeric vector specifying the limits for the x-axis.}
  \item{ylim}{a numeric vector specifying the limits for the y-axis.}
  \item{xlab}{a string giving the label for the x-axis.}
  \item{ylab}{a string giving the label for the y-axis.}
  \item{restore.par}{a logical value indicating whether the previous \code{par} settings must be restored. If you need to add points, lines, etc. to a control chart set this to \code{FALSE}.}
  \item{\dots}{additional arguments to be passed to the generic \code{\link{points}} function.} 
}
%\details{}
%\value{}
\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM. \cr
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons. \cr
Ryan, T. P. (2000), \emph{Statistical Methods for Quality Improvement}, 2nd ed. New York: John Wiley & Sons, Inc.
}
\author{Luca Scrucca}
%\note{}

\seealso{\code{\link{mqcc}}, \code{\link{stats.T2}}, \code{\link{stats.T2.single}}}
\examples{
# See examples in help(mqcc)
}

\keyword{htest}
\keyword{hplot}
\keyword{multivariate}
