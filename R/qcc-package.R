# TODO: replace @import with @importFrom for best practice
# TODO: Review that the @importFrom arguments are actually still in use
# TODO: When documenting a method, always include a link back to the generic using [generic_name()] so the reader can easily find the full documentation and other methods.

#' Quality Control Charts
#'
#' Shewhart quality control charts for continuous, attribute and count data.
#' Cusum and EWMA charts. Operating characteristic curves. Process capability
#' analysis. Pareto chart and cause-and-effect chart. Multivariate control
#' charts.
#'
#' See [vignette and documentation](../doc/index.html) accompanying the
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
#' @seealso [qcc()], [mqcc()], [cusum()],
#' [ewma()], [ocCurves()], [processCapability()],
#' [paretoChart()], [causeEffectDiagram()].
#' @references Scrucca, L. (2004). qcc: an R package for quality control
#' charting and statistical process control. *R News* 4/1, 11-17.
#' @keywords package
NULL


# TODO: Add this to the indivdual depracated functions in R/depracated.R

#' Deprecated Functions in qcc package
#'
#' These functions are provided for compatibility with older versions of the
#' package qcc, and they will likely be removed in the future.
#' 
#' @aliases qcc.groups pareto.chart process.capability oc.curves cause.and.effect qcc.overdispersion.test
#' @export qcc.groups
#' @export pareto.chart
#' @export process.capability
#' @export oc.curves
#' @export cause.and.effect
#' @export qcc.overdispersion.test
#' @param ... all arguments are passed down.
#' @seealso [deprecated()]
NULL



#' Internal 'qcc' functions
#'
#' Internal functions for package qcc.
#'
#' These functions are not intended to be called by the user.
#' 
#' @aliases .printShortMatrix print.qccplot .qcc.options qcc.c4 qccStartupMessage
#' @author Luca Scrucca
#' @keywords internal package
NULL
