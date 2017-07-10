\name{qcc}
\alias{qcc}
\alias{print.qcc}
\alias{summary.qcc}
\alias{plot.qcc}

\title{Quality Control Charts}

\description{Create an object of class 'qcc' to perform statistical quality control. This object may then be used to plot Shewhart charts, drawing OC curves, computes capability indices, and more.}

\usage{
qcc(data, type, sizes, center, std.dev, limits,
    data.name, labels, newdata, newsizes, newdata.name,
    newlabels, nsigmas = 3, confidence.level,
    rules = shewhart.rules, plot = TRUE, \dots)

\method{print}{qcc}(x, \dots)

\method{summary}{qcc}(object, digits = getOption("digits"), \dots)

\method{plot}{qcc}(x, add.stats = TRUE, chart.all = TRUE, 
     label.limits = c("LCL ", "UCL"), title, xlab, ylab, ylim, 
     axes.las = 0, digits = getOption("digits"),
     restore.par = TRUE, \dots)
}
\arguments{

  \item{data}{a data frame, a matrix or a vector containing observed data for
  the variable to chart. Each row of a data frame or a matrix, and each
  value of a vector, refers to a sample or ''rationale group''.}

  \item{type}{a character string specifying the group statistics to compute. \cr
  Available methods are:

\tabular{lll}{
              \tab Statistic charted   \tab Chart description \cr
\code{"xbar"} \tab mean                \tab means of a continuous process variable \cr
\code{"R"}    \tab range               \tab ranges of a continuous process variable \cr
\code{"S"}    \tab standard deviation  \tab standard deviations of a continuous variable \cr
\code{"xbar.one"} \tab mean            \tab one-at-time data of a continuous process variable \cr
\code{"p"}    \tab proportion          \tab proportion of nonconforming units \cr
\code{"np"}   \tab count               \tab number of nonconforming units \cr
\code{"c"}    \tab count               \tab nonconformities per unit \cr
\code{"u"}    \tab count               \tab average nonconformities per unit \cr
\code{"g"}    \tab count               \tab number of non-events between events
}

Furthermore, a user specified type of chart, say \code{"newchart"}, can be provided. This requires the definition of \code{"stats.newchart"}, \code{"sd.newchart"}, and \code{"limits.newchart"}. As an example, see \code{\link{stats.xbar}}.
}

\item{sizes}{a value or a vector of values specifying the sample sizes associated with each group. For continuous data provided as data frame or matrix the sample sizes are obtained counting the non-\code{NA} elements of each row. For \code{"p"}, \code{"np"} and \code{"u"} charts the argument \code{sizes} is required.}

\item{center}{a value specifying the center of group statistics or the ''target'' value of the process.}

\item{std.dev}{a value or an available method specifying the within-group standard deviation(s) of the process. \cr
Several methods are available for estimating the standard deviation in case of a continuous process variable; see \code{\link{sd.xbar}}, \code{\link{sd.xbar.one}}, \code{\link{sd.R}}, \code{\link{sd.S}}.
}

\item{limits}{a two-values vector specifying control limits.}

\item{data.name}{a string specifying the name of the variable which appears on the plots. If not provided is taken from the object given as data.}

\item{labels}{a character vector of labels for each group.}

\item{newdata}{a data frame, matrix or vector, as for the \code{data} argument, providing further data to plot but not included in the computations.}

\item{newsizes}{a vector as for the \code{sizes} argument providing further data sizes to plot but not included in the computations.}

\item{newdata.name}{a string specifying the name of the variable which appears on the plots. If not provided is taken from the object given as newdata.}

\item{newlabels}{a character vector of labels for each new group defined in the argument \code{newdata}.}

\item{nsigmas}{a numeric value specifying the number of sigmas to use for computing control limits. It is ignored when the \code{confidence.level} argument is  provided.} 

\item{confidence.level}{a numeric value between 0 and 1 specifying the confidence level of the computed probability limits.} 

\item{rules}{a function of rules to apply to the chart. By default, the \code{shewhart.rules} function is used.}

\item{plot}{logical. If \code{TRUE} a Shewhart chart is plotted.}

\item{add.stats}{a logical value indicating whether statistics and other information should be printed at the bottom of the chart.}

\item{chart.all}{a logical value indicating whether both statistics for \code{data} and for \code{newdata} (if given) should be plotted.}

\item{label.limits}{a character vector specifying the labels for control limits.}

\item{title}{a string giving the label for the main title.}

\item{xlab}{a string giving the label for the x-axis.}
 
\item{ylab}{a string giving the label for the y-axis.}

\item{ylim}{a numeric vector specifying the limits for the y-axis.}

\item{axes.las}{numeric in \{0,1,2,3\} specifying the style of axis labels. See \code{help(par)}.}

\item{digits}{the number of significant digits to use.}

\item{restore.par}{a logical value indicating whether the previous \code{par} settings must be restored. If you need to add points, lines, etc. to a control chart set this to \code{FALSE}.}

\item{object}{an object of class 'qcc'.}

\item{x}{an object of class 'qcc'.}

\item{\dots}{additional arguments to be passed to the generic function.} 
}
%\details{}

\value{Returns an object of class 'qcc'.}

\note{For a nice blog post discussing the \code{qcc} package, in particular how to implement the \emph{Western Eletric Rules (WER)}, see \url{http://blog.yhathq.com/posts/quality-control-in-r.html}.}

\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM. \cr
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons. \cr
Ryan, T. P. (2000), \emph{Statistical Methods for Quality Improvement}, 2nd ed. New York: John Wiley & Sons, Inc. \cr
Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}

\author{Luca Scrucca}
%\note{ ~~further notes~~ }

\seealso{\code{\link{shewhart.rules}}, \code{\link{cusum}}, \code{\link{ewma}}, \code{\link{oc.curves}}, \code{\link{process.capability}}, \code{\link{qcc.groups}} }

\examples{
##
##  Continuous data 
##
data(pistonrings)
attach(pistonrings)
diameter <- qcc.groups(diameter, sample)

qcc(diameter[1:25,], type="xbar")
qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,])
q <- qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], plot=FALSE)
plot(q, chart.all=FALSE)
qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], nsigmas=2)
qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], confidence.level=0.99)

qcc(diameter[1:25,], type="R")
qcc(diameter[1:25,], type="R", newdata=diameter[26:40,])

qcc(diameter[1:25,], type="S")
qcc(diameter[1:25,], type="S", newdata=diameter[26:40,])

# add warning limits at 2 std. deviations
q <- qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], plot=FALSE)
(warn.limits <- limits.xbar(q$center, q$std.dev, q$sizes, 2))
plot(q, restore.par = FALSE)
abline(h = warn.limits, lty = 3, col = "chocolate")

# variable control limits
out <- c(9, 10, 30, 35, 45, 64, 65, 74, 75, 85, 99, 100)
diameter <- qcc.groups(pistonrings$diameter[-out], sample[-out])
qcc(diameter[1:25,], type="xbar")
qcc(diameter[1:25,], type="R")
qcc(diameter[1:25,], type="S")
qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,])
qcc(diameter[1:25,], type="R", newdata=diameter[26:40,])
qcc(diameter[1:25,], type="S", newdata=diameter[26:40,])

detach(pistonrings)

##
##  Attribute data 
##

data(orangejuice)
attach(orangejuice)
qcc(D[trial], sizes=size[trial], type="p")

# remove out-of-control points (see help(orangejuice) for the reasons)
inc <- setdiff(which(trial), c(15,23))
q1 <- qcc(D[inc], sizes=size[inc], type="p")
qcc(D[inc], sizes=size[inc], type="p", newdata=D[!trial], newsizes=size[!trial]) 
detach(orangejuice)

data(orangejuice2)
attach(orangejuice2)
names(D) <- sample
qcc(D[trial], sizes=size[trial], type="p")
q2 <- qcc(D[trial], sizes=size[trial], type="p", newdata=D[!trial], newsizes=size[!trial])
detach(orangejuice2)

# put on the same graph the two orange juice samples
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(1,2), mar=c(5,5,3,0))
plot(q1, title="First samples", ylim=c(0,0.5), add.stats=FALSE, restore.par=FALSE)
par("mar"=c(5,0,3,3), yaxt="n")
plot(q2, title="Second samples", add.stats=FALSE, ylim=c(0,0.5))
par(oldpar)

data(circuit)
attach(circuit)
qcc(x[trial], sizes=size[trial], type="c")
# remove out-of-control points (see help(circuit) for the reasons)
inc <- setdiff(which(trial), c(6,20))
qcc(x[inc], sizes=size[inc], type="c", labels=inc)
qcc(x[inc], sizes=size[inc], type="c", labels=inc, 
    newdata=x[!trial], newsizes=size[!trial], newlabels=which(!trial))
qcc(x[inc], sizes=size[inc], type="u", labels=inc, 
    newdata=x[!trial], newsizes=size[!trial], newlabels=which(!trial))
detach(circuit)

data(pcmanufact)
attach(pcmanufact)
qcc(x, sizes=size, type="u")
detach(pcmanufact)

data(dyedcloth)
attach(dyedcloth)
qcc(x, sizes=size, type="u")
# standardized control chart
q <- qcc(x, sizes=size, type="u", plot=FALSE)
z <- (q$statistics - q$center)/sqrt(q$center/q$size)
plot(z,  type="o", ylim=range(z,3,-3), pch=16)
abline(h=0, lty=2)
abline(h=c(-3,3), lty=2)
detach(dyedcloth)

##
##  Continuous one-at-time data 
##

# viscosity data (Montgomery, pag. 242)
x <- c(33.75, 33.05, 34, 33.81, 33.46, 34.02, 33.68, 33.27, 33.49, 33.20,
       33.62, 33.00, 33.54, 33.12, 33.84)
qcc(x, type="xbar.one")
qcc(x, type="xbar.one", std.dev = "SD")
}
\keyword{htest}
\keyword{hplot}
