\name{ocCurves}
\alias{ocCurves}
\alias{print.ocCurves}
\alias{plot.ocCurves}
\alias{ocCurves.xbar}
\alias{ocCurves.R}
\alias{ocCurves.S}
\alias{ocCurves.p}
\alias{ocCurves.c}

\title{Operating Characteristic Function}

\description{Draws the operating characteristic curves for a \code{'qcc'} object.}

\usage{
ocCurves(object, \dots)

ocCurves.xbar(object, 
              size = c(1, 5, 10, 15, 20), 
              shift = seq(0, 5, by = 0.1),
              nsigmas = object$nsigmas, \dots)

ocCurves.R(object, 
           size = c(2, 5, 10, 15, 20), 
           multiplier = seq(1, 6, by = 0.1),
           nsigmas = object$nsigmas, \dots)

ocCurves.S(object, 
           size = c(2, 5, 10, 15, 20), 
           multiplier = seq(1, 6, by = 0.1),
           nsigmas = object$nsigmas, \dots) 

ocCurves.p(object, \dots)

ocCurves.c(object, \dots)
					 
\method{print}{ocCurves}(x, digits = getOption("digits"), \dots)

\method{plot}{ocCurves}(x, what = c("beta", "ARL"), 
     title, xlab, ylab, lty, lwd, col, \dots)				 
}

\arguments{
\item{object}{an object of class \code{'qcc'}.}
\item{size}{a vector of values specifying the sample sizes for which to draw the OC curves.}
\item{shift, multiplier}{a vector of values specifying the shift or multiplier values (in units of sigma).}
\item{nsigmas}{a numeric value specifying the number of sigmas to use for computing control limits; if \code{nsigmas} is \code{NULL}, \code{object$conf} is used to set up probability limits.}
\item{x}{an object of class \code{'ocCurves'}.}
\item{digits}{the number of significant digits to use.}
\item{what}{a string specifying the quantity to plot on the y-axis. Possible values are \code{"beta"} for the probability of not detecting a shift, and \code{"ARL"} for the average run length.}
\item{title}{a character string specifying the main title. Set \code{title = NULL} to remove the title.}
\item{xlab, ylab}{a string giving the label for the x-axis and the y-axis.}
\item{lty, lwd, col}{values or vector of values controlling the line type, line width and colour of curves.}
\item{\dots}{catches further ignored arguments.}
}

\details{An operating characteristic curve graphically provides information about the probability of not detecting a shift in the process. \code{ocCurves} is a generic function which calls the proper function depending on the type of \code{'qcc'} object. Further arguments provided through \code{\dots} are passed to the specific function depending on the type of chart.

The probabilities are based on the conventional assumptions about process distributions: the normal distribution for \code{"xbar"}, \code{"R"}, and \code{"S"}, the binomial distribution for \code{"p"} and \code{"np"}, and the Poisson distribution for \code{"c"} and \code{"u"}. They are all sensitive to departures from those assumptions, but to varying degrees. The performance of the \code{"S"} chart, and especially the \code{"R"} chart, are likely to be seriously affected by longer tails.}

\value{The function returns an object of class \code{'ocCurves'} which contains a matrix or a vector of beta values (the probability of type II error) and ARL (average run length).}

\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM.

Montgomery, D.C. (2013) \emph{Introduction to Statistical Quality Control}, 7th ed. New York: John Wiley & Sons.

Ryan, T. P. (2011), \emph{Statistical Methods for Quality Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.

Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17.

Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{ ~~further notes~~ }
\seealso{\code{\link{qcc}}}
\examples{
data(pistonrings)
diameter  = qccGroups(diameter, sample, data = pistonrings)
oc  = ocCurves.xbar(qcc(diameter, type="xbar", nsigmas=3))
oc
plot(oc)

data(orangejuice)
oc  = with(orangejuice,
           ocCurves(qcc(D[trial], sizes=size[trial], type="p")))
oc
plot(oc)

data(circuit)
oc  = with(circuit,
           ocCurves(qcc(x[trial], sizes=size[trial], type="c")))
oc
plot(oc)
}
\keyword{htest}
\keyword{hplot}
