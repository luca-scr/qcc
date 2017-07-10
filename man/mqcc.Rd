\name{mqcc}
\alias{mqcc}
\alias{print.mqcc}
\alias{summary.mqcc}
\alias{plot.mqcc}
\title{Multivariate Quality Control Charts}
\description{Create an object of class 'mqcc' to perform multivariate statistical quality control.}
\usage{
mqcc(data, type = c("T2", "T2.single"), center, cov,
     limits = TRUE, pred.limits = FALSE,
     data.name, labels, newdata, newlabels, 
     confidence.level = (1 - 0.0027)^p, rules = shewhart.rules,
     plot = TRUE, \dots)

\method{print}{mqcc}(x, \dots)

\method{summary}{mqcc}(object, digits = getOption("digits"), \dots)

\method{plot}{mqcc}(x, add.stats = TRUE, chart.all = TRUE,
     label.limits = c("LCL", "UCL"), label.pred.limits = c("LPL", "UPL"),
     title, xlab, ylab, ylim, axes.las = 0, 
     digits = getOption("digits"), restore.par = TRUE, \dots)
}

\arguments{

\item{data}{For subgrouped data, a list with a data frame or a matrix for each variable to monitor. Each row of the data frame or matrix refers to a sample or ''rationale'' group. 
For individual observations, where each sample has a single observation, users can provide a list with a data frame or a matrix having a single column, or a data frame or a matrix where each rows refer to samples and columns to variables. See examples.}

\item{type}{a character string specifying the type of chart:

\tabular{ll}{
            \tab Chart description \cr
\code{"T2"} \tab Hotelling \eqn{T^2} chart for subgrouped data \cr
\code{"T2.single"} \tab Hotelling \eqn{T^2} chart for individual observations \cr
}
}

\item{center}{a vector of values to use for center of input variables.}

\item{cov}{a matrix of values to use for the covariance matrix of input variables.}

\item{limits}{a logical indicating if control limits (Phase I) must be computed (by default using \code{\link{limits.T2}} or \code{\link{limits.T2.single}}) and plotted, or a two-values vector specifying control limits.}

\item{pred.limits}{a logical indicating if prediction limits (Phase II) must be computed (by default using \code{\link{limits.T2}} or \code{\link{limits.T2.single}}) and plotted, or a two-values vector specifying prediction limits.}

\item{data.name}{a string specifying the name of the variable which appears on the plots. If not provided is taken from the object given as data.}

\item{labels}{a character vector of labels for each group.}

\item{newdata}{a data frame, matrix or vector, as for the \code{data} argument, providing further data to plot but not included in the computations.}

\item{newlabels}{a character vector of labels for each new group defined in the argument \code{newdata}.}

\item{confidence.level}{a numeric value between 0 and 1 specifying the confidence level of the computed probability limits. 
By default is set at \eqn{(1 - 0.0027)^p} where \eqn{p} is the number of variables, and \eqn{0.0027} is the probability of Type I error for a single Shewhart chart at the usual 3-sigma control level.} 

\item{rules}{a function of rules to apply to the chart. By default, the \code{shewhart.rules} function is used.}

\item{plot}{logical. If \code{TRUE} a quality chart is plotted.}

\item{add.stats}{a logical value indicating whether statistics and other information should be printed at the bottom of the chart.}

\item{chart.all}{a logical value indicating whether both statistics for \code{data} and for \code{newdata} (if given) should be plotted.}

\item{label.limits}{a character vector specifying the labels for control limits (Phase I).}

\item{label.pred.limits}{a character vector specifying the labels for prediction  control limits (Phase II).}

\item{title}{a string giving the label for the main title.}

\item{xlab}{a string giving the label for the x-axis.}
 
\item{ylab}{a string giving the label for the y-axis.}

\item{ylim}{a numeric vector specifying the limits for the y-axis.}

\item{axes.las}{numeric in \{0,1,2,3\} specifying the style of axis labels. See \code{help(par)}.}

\item{digits}{the number of significant digits to use when \code{add.stats = TRUE}.}

\item{restore.par}{a logical value indicating whether the previous \code{par} settings must be restored. If you need to add points, lines, etc. to a control chart set this to \code{FALSE}.}

\item{object}{an object of class 'mqcc'.}

\item{x}{an object of class 'mqcc'.}

\item{\dots}{additional arguments to be passed to the generic function.} 
}

%\details{}
\value{Returns an object of class 'mqcc'.}

\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM. \cr
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons. \cr
Ryan, T. P. (2000), \emph{Statistical Methods for Quality Improvement}, 2nd ed. New York: John Wiley & Sons, Inc. \cr
Scrucca, L. (2004). qcc: an R package for quality control charting and statistical process control. \emph{R News} 4/1, 11-17. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{}

\seealso{\code{\link{stats.T2}}, \code{\link{stats.T2.single}}, \code{\link{limits.T2}}, \code{\link{limits.T2.single}}, \code{\link{ellipseChart}}, \code{\link{qcc}}
}

\examples{
##
##  Subgrouped data
##

# Ryan (2000, Table 9.2) data with p = 2 variables, m = 20 samples, n = 4 sample size:
X1 <- matrix(c(72, 56, 55, 44, 97, 83, 47, 88, 57, 26, 46, 
49, 71, 71, 67, 55, 49, 72, 61, 35, 84, 87, 73, 80, 26, 89, 66, 
50, 47, 39, 27, 62, 63, 58, 69, 63, 51, 80, 74, 38, 79, 33, 22, 
54, 48, 91, 53, 84, 41, 52, 63, 78, 82, 69, 70, 72, 55, 61, 62, 
41, 49, 42, 60, 74, 58, 62, 58, 69, 46, 48, 34, 87, 55, 70, 94, 
49, 76, 59, 57, 46), ncol = 4)
X2 <- matrix(c(23, 14, 13, 9, 36, 30, 12, 31, 14, 7, 10, 
11, 22, 21, 18, 15, 13, 22, 19, 10, 30, 31, 22, 28, 10, 35, 18, 
11, 10, 11, 8, 20, 16, 19, 19, 16, 14, 28, 20, 11, 28, 8, 6, 
15, 14, 36, 14, 30, 8, 35, 19, 27, 31, 17, 18, 20, 16, 18, 16, 
13, 10, 9, 16, 25, 15, 18, 16, 19, 10, 30, 9, 31, 15, 20, 35, 
12, 26, 17, 14, 16), ncol = 4)
X <- list(X1 = X1, X2 = X2)

q <- mqcc(X, type = "T2")
summary(q)
ellipseChart(q)
ellipseChart(q, show.id = TRUE)

q <- mqcc(X, type = "T2", pred.limits = TRUE)

# Ryan (2000) discussed Xbar-charts for single variables computed adjusting the 
# confidence level of the T^2 chart:
q1 <- qcc(X1, type = "xbar", confidence.level = q$confidence.level^(1/2))
summary(q1)
q2 <- qcc(X2, type = "xbar", confidence.level = q$confidence.level^(1/2))
summary(q2)


require(MASS)
# generate new "in control" data
Xnew <- list(X1 = matrix(NA, 10, 4), X2 =  matrix(NA, 10, 4))
for(i in 1:4)
   { x <- mvrnorm(10, mu = q$center, Sigma = q$cov)
     Xnew$X1[,i] <- x[,1]
     Xnew$X2[,i] <- x[,2] 
   }
qq <- mqcc(X, type = "T2", newdata = Xnew, pred.limits = TRUE)
summary(qq)

# generate new "out of control" data
Xnew <- list(X1 = matrix(NA, 10, 4), X2 =  matrix(NA, 10, 4))
for(i in 1:4)
   { x <- mvrnorm(10, mu = 1.2*q$center, Sigma = q$cov)
     Xnew$X1[,i] <- x[,1]
     Xnew$X2[,i] <- x[,2] 
   }
qq <- mqcc(X, type = "T2", newdata = Xnew, pred.limits = TRUE)
summary(qq)

##
## Individual observations data
##

data(boiler)

q <- mqcc(boiler, type = "T2.single", confidence.level = 0.999)
summary(q)

# generate new "in control" data
boilerNew <- mvrnorm(10, mu = q$center, Sigma = q$cov)
qq <- mqcc(boiler, type = "T2.single", confidence.level = 0.999, 
           newdata = boilerNew, pred.limits = TRUE)
summary(qq)

# generate new "out of control" data
boilerNew = mvrnorm(10, mu = 1.01*q$center, Sigma = q$cov)
qq <- mqcc(boiler, type = "T2.single", confidence.level = 0.999, 
           newdata = boilerNew, pred.limits = TRUE)
summary(qq)

# provides "robust" estimates of means and covariance matrix
library(MASS)
rob <- cov.rob(boiler)
qrob <- mqcc(boiler, type = "T2.single", center = rob$center, cov = rob$cov)
summary(qrob)

}

\keyword{htest}
\keyword{hplot}
\keyword{multivariate}
