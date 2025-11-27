\name{qcc}
\alias{qcc}
\alias{print.qcc}
\alias{summary.qcc}
\alias{plot.qcc}

\title{Quality Control Charts}

\description{Create an object of class \code{'qcc'} to perform statistical quality control. This object may then be used to plot Shewhart charts, drawing OC curves, computes capability indices, and more.}

\usage{
qcc(data, 
    type = c("xbar", "R", "S", "xbar.one", "p", "np", "c", "u", "g"), 
    sizes, center, std.dev, limits,
    newdata, newsizes, 
    nsigmas = 3, confidence.level, 
    rules = c(1,4), \dots)

\method{print}{qcc}(x, digits = getOption("digits"), \dots)

\method{plot}{qcc}(x,
     xtime = NULL,
     add.stats = qcc.options("add.stats"), 
     chart.all = qcc.options("chart.all"), 
     fill = qcc.options("fill"),
     label.center = "CL",
     label.limits = c("LCL ", "UCL"),
     title, xlab, ylab, xlim, ylim,
     digits = getOption("digits"), \dots)
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
\code{"g"}    \tab count               \tab number of non-events between events \cr
}
Furthermore, a user specified type of chart, say \code{"newchart"}, can be provided. This requires the definition of \code{"stats.newchart"}, \code{"sd.newchart"}, and \code{"limits.newchart"}. As an example, see \code{\link{stats.xbar}}.
}

\item{sizes}{a value or a vector of values specifying the sample sizes associated with each group. For continuous data provided as data frame or matrix the sample sizes are obtained counting the non-\code{NA} elements of each row. For \code{"p"}, \code{"np"} and \code{"u"} charts the argument \code{sizes} is required.}

\item{center}{a value specifying the center of group statistics or the ''target'' value of the process.}

\item{std.dev}{a value or an available method specifying the within-group standard deviation(s) of the process. \cr
Several methods are available for estimating the standard deviation in case of a continuous process variable; see \code{\link{sd.xbar}}, \code{\link{sd.xbar.one}}, \code{\link{sd.R}}, \code{\link{sd.S}}.
}

\item{limits}{a two-values vector specifying control limits.}

\item{newdata}{a data frame, matrix or vector, as for the \code{data} argument, providing further data to plot but not included in the computations.}

\item{newsizes}{a vector as for the \code{sizes} argument providing further data sizes to plot but not included in the computations.}

\item{nsigmas}{a numeric value specifying the number of sigmas to use for computing control limits. It is ignored when the \code{confidence.level} argument is  provided.} 

\item{confidence.level}{a numeric value between 0 and 1 specifying the confidence level of the computed probability limits.} 

\item{rules}{a value or a vector of values specifying the rules to apply to the chart. See \code{\link{qccRules}} for possible values and their meaning.}

\item{xtime}{a vector of date-time values as returned by \code{\link{Sys.time}} and \code{\link{Sys.Date}}. If provided it is used for x-axis so it must be of the same length as the statistic charted.}

\item{add.stats}{a logical value indicating whether statistics and other information should be printed at the bottom of the chart.}

\item{chart.all}{a logical value indicating whether both statistics for \code{data} and for \code{newdata} (if given) should be plotted.}

\item{fill}{a logical value specifying if the in-control area should be filled with the color specified in \code{qcc.options("zones")$fill}.}

\item{label.center}{a character  specifying the label for center line.}

\item{label.limits}{a character vector specifying the labels for control limits.}

\item{title}{a character string specifying the main title. Set \code{title = NULL} to remove the title.}

\item{xlab, ylab}{a string giving the label for the x-axis and the y-axis.}

\item{xlim, ylim}{a numeric vector specifying the limits for the x-axis and the y-axis.}

\item{digits}{the number of significant digits to use.}

\item{x}{an object of class \code{'qcc'}.}

\item{\dots}{additional arguments to be passed to the generic function.} 
}
%\details{}

\value{Returns an object of class \code{'qcc'}.}

\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM.

Montgomery, D.C. (2013) \emph{Introduction to Statistical Quality Control}, 7th ed. New York: John Wiley & Sons.

Ryan, T. P. (2011), \emph{Statistical Methods for Quality Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.

Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17. 

Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}

\author{Luca Scrucca}
%\note{ ~~further notes~~ }

\seealso{\code{\link{qccRules}}, \code{\link{cusum}}, \code{\link{ewma}}, \code{\link{ocCurves}}, \code{\link{processCapability}}, \code{\link{qccGroups}} }

\examples{
##
##  Continuous data 
##
data(pistonrings)
diameter  = qccGroups(data = pistonrings, diameter, sample)

(q  = qcc(diameter[1:25,], type="xbar"))
plot(q)

(q  = qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,]))
plot(q)

q  = qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,])
plot(q, chart.all=FALSE)

plot(qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], nsigmas=2))

plot(qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], confidence.level=0.99))

q  = qcc(diameter[1:25,], type="R")
q
plot(q)

plot(qcc(diameter[1:25,], type="R", newdata=diameter[26:40,]))

plot(qcc(diameter[1:25,], type="S"))

plot(qcc(diameter[1:25,], type="S", newdata=diameter[26:40,]))

plot(qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,], rules = 1:4))

# variable control limits
out  = c(9, 10, 30, 35, 45, 64, 65, 74, 75, 85, 99, 100)
diameter  = qccGroups(data = pistonrings[-out,], diameter, sample)
plot(qcc(diameter[1:25,], type="xbar"))
plot(qcc(diameter[1:25,], type="R"))
plot(qcc(diameter[1:25,], type="S"))
plot(qcc(diameter[1:25,], type="xbar", newdata=diameter[26:40,]))
plot(qcc(diameter[1:25,], type="R", newdata=diameter[26:40,]))
plot(qcc(diameter[1:25,], type="S", newdata=diameter[26:40,]))

# to customize a Shewhart chart use
q = qcc(diameter[1:25, ], type = "xbar")
plot = plot(q, ylim = c(73.9, 74.1)) 
# returned object is of class "patchwork", then add geom_* layer to the first element
plot[[1]] <- plot[[1]] + 
  geom_hline(yintercept = c(73.95, 74.05), lty = 2) 
plot

##
##  Attribute data 
##

data(orangejuice)

q  = with(orangejuice, 
          qcc(D[trial], sizes=size[trial], type="p"))
q
plot(q)

# remove out-of-control points (see help(orangejuice) for the reasons)
outofctrl  = c(15,23)
q1  = with(orangejuice[-outofctrl,], 
           qcc(D[trial], sizes=size[trial], type="p"))
plot(q1)
q1  = with(orangejuice[-outofctrl,], 
           qcc(D[trial], sizes=size[trial], type="p",
               newdata=D[!trial], newsizes=size[!trial]))
plot(q1)

data(orangejuice2)
q2  = with(orangejuice2, 
           qcc(D[trial], sizes=size[trial], type="p"))
plot(q2)
q2  = with(orangejuice2, 
           qcc(D[trial], sizes=size[trial], type="p", 
               newdata=D[!trial], newsizes=size[!trial]))
plot(q2)

data(circuit)
plot(with(circuit, qcc(x[trial], sizes=size[trial], type="c")))

# remove out-of-control points (see help(circuit) for the reasons)
outofctrl  = c(15,23)
q1  = with(orangejuice[-outofctrl,], 
           qcc(D[trial], sizes=size[trial], type="p"))
plot(q1)
q1  = with(orangejuice[-outofctrl,], 
           qcc(D[trial], sizes=size[trial], type="p",
               newdata=D[!trial], newsizes=size[!trial]))
plot(q1)

outofctrl  = c(6,20)
q1  = with(circuit[-outofctrl,], 
           qcc(x[trial], sizes=size[trial], type="c"))
plot(q1)
q1  = with(circuit[-outofctrl,], 
           qcc(x[trial], sizes=size[trial], type="c", 
               newdata = x[!trial], newsizes = size[!trial]))
plot(q1)
q1  = with(circuit[-outofctrl,], 
           qcc(x[trial], sizes=size[trial], type="u", 
           newdata = x[!trial], newsizes = size[!trial]))
plot(q1)

data(pcmanufact)
q1  = with(pcmanufact, qcc(x, sizes=size, type="u"))
q1
plot(q1)

data(dyedcloth)
# variable control limits
plot(with(dyedcloth, qcc(x, sizes=size, type="u")))
# standardized control chart
q  = with(dyedcloth, qcc(x, sizes=size, type="u"))
z  = (q$statistics - q$center)/sqrt(q$center/q$size)
plot(qcc(z, sizes = 1, type = "u", center = 0, std.dev = 1, limits = c(-3,3)),
     title = "Standardized u chart")
    
##
##  Continuous one-at-time data 
##

data(viscosity)
q  = with(viscosity, 
          qcc(viscosity[trial], type = "xbar.one"))
q
plot(q)
# batch 4 is out-of-control because of a process temperature controller
# failure; remove it and recompute
viscosity  = viscosity[-4,]
plot(with(viscosity, 
          qcc(viscosity[trial], type = "xbar.one", newdata = viscosity[!trial])))
}
\keyword{htest}
\keyword{hplot}


