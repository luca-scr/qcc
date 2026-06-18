# TODO: replace @import with @importFrom for best practice
# TODO: Review that the @importFrom arguments are actually still in use


#' Quality Control Charts
#' 
#' Shewhart quality control charts for continuous, attribute and count data.
#' Cusum and EWMA charts. Operating characteristic curves. Process capability
#' analysis. Pareto chart and cause-and-effect chart. Multivariate control
#' charts.
#' 
#' See \href{../doc/index.html}{vignette and documentation} accompanying the
#' package.
#' 
#' @name qcc-package
#' @docType package
#' @import stats utils ggplot2 patchwork
#' @importFrom graphics strheight strwidth hist abline axis box contour lines mtext par points polygon rect text
#' @importFrom grDevices gray adjustcolor palette extendrange nclass.FD nclass.Sturges
#' @importFrom MASS mvrnorm
#' @importFrom scales label_percent
#' @importFrom cli rule
#' @importFrom crayon bold
#' @author Luca Scrucca
#' @seealso \code{\link{qcc}}, \code{\link{mqcc}}, \code{\link{cusum}},
#' \code{\link{ewma}}, \code{\link{ocCurves}}, \code{\link{processCapability}},
#' \code{\link{paretoChart}}, \code{\link{causeEffectDiagram}}.
#' @references Scrucca, L. (2004). qcc: an R package for quality control
#' charting and statistical process control. \emph{R News} 4/1, 11-17.
#' @keywords package
NULL



#' Deprecated Functions in \pkg{qcc} package
#' 
#' These functions are provided for compatibility with older versions of the
#' package \pkg{qcc}, and they will likely be removed in the future.
#' 
#' 
#' @aliases qcc.groups pareto.chart process.capability oc.curves
#' cause.and.effect qcc.overdispersion.test
#' @param \dots all arguments are passed down.
#' @seealso \code{\link{deprecated}}
NULL



#' Internal 'qcc' functions
#' 
#' Internal functions for package \pkg{qcc}.
#' 
#' These functions are not intended to be called by the user.
#' 
#' @aliases .printShortMatrix print.qccplot .qcc.options qcc.c4
#' qccStartupMessage
#' @author Luca Scrucca
#' @keywords internal package
NULL


