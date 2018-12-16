\name{qcc}
\alias{qcc}
\alias{print.qcc}
\alias{summary.qcc}
\alias{plot.qcc}

\title{Quality Control Charts}

\description{Create an object of class 'qcc' to perform statistical quality control. This object may then be used to plot Shewhart charts, drawing OC curves, computes capability indices, and more.}

\usage{
qcc(data, 
    type = c("xbar", "R", "S", "xbar.one", "p", "np", "c", "u", "g"), 
    sizes, center, std.dev, limits,
    data.name, labels, newdata, newsizes, newdata.name,
    newlabels, nsigmas = 3, confidence.level,
    rules = c(1,4), plot = TRUE, \dots)

\method{print}{qcc}(x, digits = getOption("digits"), \dots)

\method{plot}{qcc}(x,
	   add.stats = qcc.options("add.stats"), 
     chart.all = qcc.options("chart.all"), 
     fill = qcc.options("fill"),
     label.limits = c("LCL ", "UCL"), 
		 title, xlab, ylab, ylim, axes.las = 0, 
     digits = getOption("digits"),
     restore.par = TRUE, \dots)
}
\arguments{

  \item{data}{a data frame, a matrix or a vector containing observed data for
  the variable to chart. Each row of a data frame or a matrix, and each
  value of a vector, refers to a sample or ''rationale group''.}

  \item{type}{a character string specifying the group statistics to compute. \cr
  Available methods are:
%\tabular{lll}{
%              \tab Statistic charted   \tab Chart description \cr
%\code{"xbar"} \tab mean                \tab means of a continuous process variable \cr
%\code{"R"}    \tab range               \tab ranges of a continuous process variable \cr
%\code{"S"}    \tab standard deviation  \tab standard deviations of a continuous variable \cr
%\code{"xbar.one"} \tab mean            \tab one-at-time data of a continuous process variable \cr
%\code{"p"}    \tab proportion          \tab proportion of nonconforming units \cr
%\code{"np"}   \tab count               \tab number of nonconforming units \cr
%\code{"c"}    \tab count               \tab nonconformities per unit \cr
%\code{"u"}    \tab count               \tab average nonconformities per unit %\cr
%\code{"g"}    \tab count               \tab number of non-events between events \cr
%}

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

\item{rules}{a value or a vector of values specifying the rules to apply to the chart. See \code{\link{qccRules}} for possible values and their meaning.}

\item{plot}{logical. If \code{TRUE} a Shewhart chart is plotted.}

\item{add.stats}{a logical value indicating whether statistics and other information should be printed at the bottom of the chart.}

\item{chart.all}{a logical value indicating whether both statistics for \code{data} and for \code{newdata} (if given) should be plotted.}

\item{fill}{a logical value specifying if the in-control area should be filled with the color specified in 
\code{qcc.options("zones")$fill}.}

\item{label.limits}{a character vector specifying the labels for control limits.}

\item{title}{a string giving the label for the main title.}

\item{xlab}{a string giving the label for the x-axis.}
 
\item{ylab}{a string giving the label for the y-axis.}

\item{ylim}{a numeric vector specifying the limits for the y-axis.}

\item{axes.las}{numeric in \{0,1,2,3\} specifying the style of axis labels. See \code{help(par)}.}

\item{digits}{the number of significant digits to use.}

\item{restore.par}{a logical value indicating whether the previous \code{par} settings must be restored. If you need to add points, lines, etc. to a control chart set this to \code{FALSE}.}

\item{x}{an object of class 'qcc'.}

\item{\dots}{additional arguments to be passed to the generic function.} 
}
%\details{}

\value{Returns an object of class 'qcc'.}

\note{For a nice blog post discussing the \code{qcc} package, in particular how to implement the \emph{Western Eletric Rules (WER)}, see \url{http://blog.yhat.com/posts/quality-control-in-r.html}.}

\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM.

Montgomery, D.C. (2009) \emph{Introduction to Statistical Quality Control}, 6th ed. New York: John Wiley & Sons.

Ryan, T. P. (2011), \emph{Statistical Methods for Quality Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.

Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17. 

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
diameter <- qcc.groups(diameter, sample, data = pistonrings)

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
diameter <- qcc.groups(diameter, sample, data = pistonrings[-out,])
qcc(diameter[1:25,], type="xbar")
qcc(diameter[1:25,], type="R")
qcc(diameter[1:25,], type="S")
qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,])
qcc(diameter[1:25,], type="R", newdata=diameter[26:40,])
qcc(diameter[1:25,], type="S", newdata=diameter[26:40,])

##
##  Attribute data 
##

data(orangejuice)

with(orangejuice, 
     qcc(D[trial], sizes=size[trial], type="p"))

# remove out-of-control points (see help(orangejuice) for the reasons)
outofctrl <- c(15,23)
q1 <- with(orangejuice[-outofctrl,], 
           qcc(D[trial], sizes=size[trial], type="p"))
with(orangejuice[-outofctrl,], 
     qcc(D[trial], sizes=size[trial], type="p",
         newdata=D[!trial], newsizes=size[!trial]))

data(orangejuice2)
with(orangejuice2, 
     qcc(D[trial], sizes=size[trial], type="p", labels = sample[trial]))
q2 <- with(orangejuice2, 
           qcc(D[trial], sizes=size[trial], type="p", labels = sample[trial],
               newdata=D[!trial], newsizes=size[!trial], newlabels = sample[!trial]))

data(circuit)
with(circuit, qcc(x[trial], sizes=size[trial], type="c"))

# remove out-of-control points (see help(circuit) for the reasons)
outofctrl <- c(15,23)
q1 <- with(orangejuice[-outofctrl,], 
           qcc(D[trial], sizes=size[trial], type="p"))
with(orangejuice[-outofctrl,], 
     qcc(D[trial], sizes=size[trial], type="p",
         newdata=D[!trial], newsizes=size[!trial]))

outofctrl <- c(6,20)
with(circuit[-outofctrl,], 
     qcc(x[trial], sizes=size[trial], type="c", labels = sample[trial]))
with(circuit[-outofctrl,], 
     qcc(x[trial], sizes=size[trial], type="c", labels = sample[trial],
         newdata = x[!trial], newsizes = size[!trial], newlabels = sample[!trial]))
with(circuit[-outofctrl,], 
     qcc(x[trial], sizes=size[trial], type="u", labels = sample[trial],
         newdata = x[!trial], newsizes = size[!trial], newlabels = sample[!trial]))

data(pcmanufact)
with(pcmanufact, qcc(x, sizes=size, type="u"))

data(dyedcloth)
# variable control limits
with(dyedcloth, qcc(x, sizes=size, type="u"))
# standardized control chart
q <- with(dyedcloth, qcc(x, sizes=size, type="u", plot=FALSE))
z <- (q$statistics - q$center)/sqrt(q$center/q$size)
qcc(z, sizes = 1, type = "u", title = "Standardized u chart",
    center = 0, std.dev = 1, limits = c(-3,3))
    
##
##  Continuous one-at-time data 
##

data(viscosity)
with(viscosity, 
     qcc(viscosity[trial], type = "xbar.one", labels = batch[trial]))
# batch 4 is out-of-control because of a process temperature controller
# failure; remove it and recompute
viscosity <- viscosity[-4,]
with(viscosity, 
     qcc(viscosity[trial], type = "xbar.one", labels = batch[trial],
         newdata = viscosity[!trial], newlabels = batch[!trial]))
}
\keyword{htest}
\keyword{hplot}
