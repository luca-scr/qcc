\name{processCapability}
\alias{processCapability}
\alias{print.processCapability}
\alias{summary.processCapability}
\alias{plot.processCapability}

\title{Process capability analysis}

\description{Computes process capability indices for a \code{'qcc'} object of type \code{"xbar"} and plot the histogram.}

\usage{
processCapability(object, spec.limits, target, 
                  std.dev, nsigmas, 
                  confidence.level = 0.95,  \dots)
                   
\method{print}{processCapability}(x, digits = getOption("digits"), \dots)

\method{plot}{processCapability}(x, 
     add.stats = qcc.options("add.stats"), 
     breaks = nclass.hist, 
     fill = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.5),
     color = "white", title, xlab, 
     digits = getOption("digits"), \dots)
}
\arguments{
\item{object}{a \code{'qcc'} object of type \code{"xbar"}}
\item{spec.limits}{a two-values vector specifying the lower and upper specification limits. For one-sided specification limits, the value of the missing limit must be set to \code{NA}.}
\item{target}{a value specifying the target of the process. If missing the value from the \code{'qcc'} object is used if not \code{NULL}, otherwise the target is set at the middle value between specification limits.}
\item{std.dev}{a value specifying the within-group standard deviation. If not provided is taken from the \code{'qcc'} object.}
\item{nsigmas}{a numeric value specifying the number of sigmas to use. If not provided is taken from the \code{'qcc'} object.}
\item{confidence.level}{a numeric value between 0 and 1 specifying the level to use for computing confidence intervals.}
\item{x}{an object of class \code{'processCapability'}.}
\item{add.stats}{a logical value indicating whether statistics and capability indices should be added at the bottom of the chart.}
\item{breaks}{a value or a function used to select the number of bins in a histogram. See the help for \code{\link{nclass.scott}} for more details.}
\item{fill, color}{values specifying the colour of the filled area and the border used for drawing the histogram.}
\item{title}{a character string specifying the plot title. Set \code{title = NULL} to remove the title.}
\item{xlab}{a character string specifying the label for the x-axis.}
\item{digits}{the number of significant digits to use.}
\item{\dots}{catches further ignored arguments.}
}
\details{
This function calculates confidence limits for \eqn{C_p}{C_p} using the method described by Chou et al. (1990). 
Approximate confidence limits for \eqn{C_{pl}}{C_pl}, \eqn{C_{pu}}{C_pu} and  \eqn{C_{pk}}{C_pk} are computed using the method in Bissell (1990).  
Confidence limits for \eqn{C_{pm}}{C_pm} are based on the method of Boyles (1991); this method is approximate and it assumes that the target is midway between the specification limits.
}
\value{Invisibly returns a list with components:
\item{nobs}{number of observations}
\item{center}{center}
\item{std.dev}{standard deviation}
\item{target}{target}
\item{spec.limits}{a vector of values giving the lower specification limit (LSL) and the upper specification limit (USL)}
\item{indices}{a matrix of capability indices (\eqn{C_p}{C_p}, \eqn{C_{pl}}{C_pl}, \eqn{C_{pu}}{C_pu}, \eqn{C_{pk}}{C_pk}, \eqn{C_{pm}}{C_pm}) and the corresponding confidence limits.}
\item{exp}{a vector of values giving the expected fraction, based on a normal approximation, of the observations less than LSL and greater than USL.}
\item{obs}{a vector of values giving the fraction of observations less than LSL and greater than USL.}
}
\references{
Bissell, A.F. (1990) \emph{How reliable is your capability index?}, Applied Statistics, 39, 331-340.

Boyles, R.A. (1991) \emph{The Taguchi capability index}, Journal of Quality Technology, 23, 107-126.

Chou, Y., Owen D.B. and Borrego S.A. (1990) \emph{Lower Confidence Limits on Process Capability Indices}, Journal of Quality Technology, 22, 223-229.

Montgomery, D.C. (2013) \emph{Introduction to Statistical Quality Control}, 7th ed. New York: John Wiley & Sons.

Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
\seealso{\code{\link{qcc}}}
\examples{
data(pistonrings)
diameter  = qccGroups(data = pistonrings, diameter, sample)
q  = qcc(diameter[1:25,], type="xbar", nsigmas=3)
pc  = processCapability(q, spec.limits=c(73.95,74.05))
pc
plot(pc)
plot(processCapability(q, spec.limits=c(73.95,74.05), target=74.02))
plot(processCapability(q, spec.limits=c(73.99,74.01)))
plot(processCapability(q, spec.limits = c(73.99, 74.1)))
}
\keyword{htest}
\keyword{hplot}

