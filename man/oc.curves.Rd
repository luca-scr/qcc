\name{ocCurves}
\alias{ocCurves}
\alias{ocCurves.xbar}
\alias{ocCurves.R}
\alias{ocCurves.S}
\alias{ocCurves.p}
\alias{ocCurves.c}

\title{Operating Characteristic Function}

\description{Draws the operating characteristic curves for a \code{'qcc'} object.}

\usage{
ocCurves(object, \dots)

ocCurves.xbar(object, n, c = seq(0, 5, length=101), nsigmas = object$nsigmas, 
              lty = rep(1, length(n)), lwd = rep(2, length(n)), 
              col = blues.colors(length(n)), 
              identify=FALSE, restore.par=TRUE)

ocCurves.R(object, n, c = seq(1, 6, length=101), nsigmas = object$nsigmas, 
           lty = rep(1, length(n)), lwd = rep(2, length(n)), 
           col = blues.colors(length(n)), 
           identify = FALSE, restore.par=TRUE)

ocCurves.S(object, n, c = seq(1, 6, length=101), nsigmas = object$nsigmas, 
           lty = rep(1, length(n)), lwd = rep(2, length(n)), 
           col = blues.colors(length(n)), 
           identify = FALSE, restore.par=TRUE)

ocCurves.p(object, nsigmas = object$nsigmas, 
           lty = 1, lwd = 2, col = blues.colors(1),
           identify = FALSE, restore.par=TRUE)

ocCurves.c(object, nsigmas = object$nsigmas, 
           lty = 1, lwd = 2, col = blues.colors(1),
           identify = FALSE, restore.par=TRUE)
}

\arguments{
\item{object}{an object of class \code{'qcc'}.}
\item{identify}{logical specifying whether to interactively identify points on the plot (see help for \code{\link{identify}}).}
\item{n}{a vector of values specifying the sample sizes for which to draw the OC curves.}
\item{c}{a vector of values specifying the multipliers for sigma in case of continuous variable.}
\item{nsigmas}{a numeric value specifying the number of sigmas to use for
computing control limits; if \code{nsigmas} is \code{NULL},
\code{object$conf} is used to set up probability limits; \code{nsigmas}
is ignored for types "p" and "c".}
\item{lty, lwd, col}{values or vector of values controlling the line type, line width and colour of curves.}
\item{restore.par}{a logical value indicating whether the previous \code{par} settings must be restored. If you need to add points, lines, etc. to a chart set this to \code{FALSE}.}
\item{\dots}{additional arguments to be passed to the generic function.}
}
\details{An operating characteristic curve graphically provides information about the probability of not detecting a shift in the process. \code{ocCurves} is a generic function which calls the proper function depending on the type of \code{'qcc'} object. Further arguments provided through \code{\dots} are passed to the specific function depending on the type of chart.

The probabilities are based on the conventional assumptions about process distributions: the normal distribution for \code{"xbar"}, \code{"R"}, and \code{"S"}, the binomial distribution for \code{"p"} and \code{"np"}, and the Poisson distribution for \code{"c"} and \code{"u"}. They are all sensitive to departures from those assumptions, but to varying degrees. The performance of the \code{"S"} chart, and especially the \code{"R"} chart, are likely to be seriously affected by longer tails.}

\value{The function invisibly returns a matrix or a vector of beta values, the probability of type II error.}

\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM.

Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons.

Ryan, T. P. (2011), \emph{Statistical Methods for Quality Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.

Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17.

Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{ ~~further notes~~ }
\seealso{\code{\link{qcc}}}
\examples{
data(pistonrings)
diameter <- qccGroups(diameter, sample, data = pistonrings)
beta <- ocCurves.xbar(qcc(diameter, type="xbar", nsigmas=3, plot=FALSE))
print(round(beta, digits=4))
# or to identify points on the plot use
\dontrun{ocCurves.xbar(qcc(diameter, type="xbar", nsigmas=3, plot=FALSE), identify=TRUE)}

data(orangejuice)
beta <- with(orangejuice,
             ocCurves(qcc(D[trial], sizes=size[trial], type="p", plot=FALSE)))
print(round(beta, digits=4))

data(circuit)
q <- with(circuit,
          qcc(x[trial], sizes=size[trial], type="c", plot=FALSE))
beta <- ocCurves(q)
print(round(beta, digits=4))
}
\keyword{htest}
\keyword{hplot}
