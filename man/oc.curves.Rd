\name{oc.curves}
\alias{oc.curves}
\alias{oc.curves.xbar}
\alias{oc.curves.R}
\alias{oc.curves.S}
\alias{oc.curves.p}
\alias{oc.curves.c}
\title{Operating Characteristic Function}
\description{Draws the operating characteristic curves for a 'qcc' object.}
\usage{
oc.curves(object, \dots)

oc.curves.xbar(object, n, c = seq(0, 5, length=101), 
               nsigmas = object$nsigmas, identify=FALSE, restore.par=TRUE)

oc.curves.R(object, n, c = seq(1, 6, length=101),
            nsigmas = object$nsigmas, identify = FALSE, restore.par=TRUE)

oc.curves.S(object, n, c = seq(1, 6, length=101),
            nsigmas = object$nsigmas, identify = FALSE, restore.par=TRUE)

oc.curves.p(object, nsigmas = object$nsigmas, identify = FALSE, restore.par=TRUE)

oc.curves.c(object, nsigmas = object$nsigmas, identify = FALSE, restore.par=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{an object of class 'qcc'.}
\item{identify}{logical specifying whether to interactively identify points on the plot (see help for \code{\link{identify}}).}
\item{n}{a vector of values specifying the sample sizes for which to draw the OC curves.}
\item{c}{a vector of values specifying the multipliers for sigma in case of continuous variable.}
\item{nsigmas}{a numeric value specifying the number of sigmas to use for
computing control limits; if \code{nsigmas} is \code{NULL},
\code{object$conf} is used to set up probability limits; \code{nsigmas}
is ignored for types "p" and "c".}
\item{restore.par}{a logical value indicating whether the previous \code{par} settings must be restored. If you need to add points, lines, etc. to a chart set this to \code{FALSE}.}
\item{\dots}{additional arguments to be passed to the generic function.}
}
\details{An operating characteristic curve graphically provides information about the probability of not detecting a shift in the process. \code{oc.curves} is a generic function which calls the proper function depending on the type of 'qcc' object. Further arguments provided through \dots are passed to the specific function depending on the type of chart.

The probabilities are based on the conventional assumptions about
process distributions: the normal distribution for "xbar" , "R", and "S",
the binomial distribution for "p" and "np", and the Poisson distribution
for "c" and "u". They are all sensitive to departures from those assumptions,
but to varying degrees. The performance of the "S" chart, and especially
the "R" chart,
are likely to be seriously affected by longer tails.}
\value{The function invisibly returns a matrix or a vector of beta values, the probability of type II error.}
\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM. \cr
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons. \cr
Ryan, T. P. (2000), \emph{Statistical Methods for Quality Improvement}, 2nd ed. New York: John Wiley & Sons, Inc. \cr
Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{ ~~further notes~~ }
\seealso{\code{\link{qcc}}}
\examples{
data(pistonrings)
attach(pistonrings)
diameter <- qcc.groups(diameter, sample)
beta <- oc.curves.xbar(qcc(diameter, type="xbar", nsigmas=3, plot=FALSE))
print(round(beta, digits=4))
# or to identify points on the plot use
\dontrun{oc.curves.xbar(qcc(diameter, type="xbar", nsigmas=3, plot=FALSE), identify=TRUE)}
detach(pistonrings)

data(orangejuice)
attach(orangejuice)
beta <- oc.curves(qcc(D[trial], sizes=size[trial], type="p", plot=FALSE))
print(round(beta, digits=4))
# or to identify points on the plot use
\dontrun{oc.curves(qcc(D[trial], sizes=size[trial], type="p", plot=FALSE), identify=TRUE)}
detach(orangejuice)

data(circuit)
attach(circuit)
q <- qcc(x[trial], sizes=size[trial], type="c", plot=FALSE)
beta <- oc.curves(q)
print(round(beta, digits=4))
# or to identify points on the plot use
\dontrun{oc.curves(qcc(x[trial], sizes=size[trial], type="c", plot=FALSE), identify=TRUE)}
detach(circuit)
}
\keyword{htest}
\keyword{hplot}
