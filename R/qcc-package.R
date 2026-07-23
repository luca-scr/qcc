
#----------------------------------------------------------------------------#
#                                                                            #
#                     QUALITY CONTROL CHARTS IN R                            #
#                                                                            #
#  An R package for statistical in-line quality control.                     #
#                                                                            #
#  Written by: Luca Scrucca                                                  #
#              Department of Economics                                       #
#              University of Perugia, ITALY                                  #
#              luca.scrucca@unipg.it                                         #
#                                                                            #
#----------------------------------------------------------------------------#

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
#' @author Luca Scrucca
#' @seealso [qcc()], [mqcc()], [cusum()],
#' [ewma()], [ocCurves()], [processCapability()],
#' [paretoChart()], [causeEffectDiagram()].
#' @references `r refs("scrucca_2004")`
#' @keywords package internal
"_PACKAGE"

## usethis namespace: start
#' @import stats ggplot2 patchwork
#' @importFrom utils packageVersion
#' @importFrom graphics strheight strwidth hist abline axis box contour lines mtext par points polygon rect text
#' @importFrom grDevices gray adjustcolor palette extendrange nclass.FD nclass.Sturges
#' @importFrom MASS mvrnorm
#' @importFrom scales label_percent
#' @importFrom cli rule
#' @importFrom crayon bold
## usethis namespace: end
NULL

