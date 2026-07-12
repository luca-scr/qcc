#' Deprecated Functions in qcc package
#'
#' These functions are provided for compatibility with older versions of the
#' package qcc, and they will likely be removed in the future.
#'
#' @param ... all arguments are passed down.
#' @seealso [base::.Deprecated()]
#' @aliases qcc-deprecated
#' @rdname qcc-deprecated
#' @export
pareto.chart <- function(...) 
{
  .Deprecated("paretoChart", package = "qcc")
  paretoChart(...)
}

#' @rdname qcc-deprecated
#' @export
process.capability <- function(...) 
{
  .Deprecated("processCapability", package = "qcc")
  processCapability(...)
}

#' @rdname qcc-deprecated
#' @export
oc.curves <- function(...) 
{
  .Deprecated("ocCurves", package = "qcc")
  ocCurves(...)
}

#' @rdname qcc-deprecated
#' @export
qcc.overdispersion.test <- function(...) 
{
  .Deprecated("qccOverdispersionTest", package = "qcc")
  qccOverdispersionTest(...)
}

#' @rdname qcc-deprecated
#' @export
cause.and.effect <- function(...) 
{
  .Deprecated("causeEffectDiagram", package = "qcc")
  causeEffectDiagram(...)
}

#' @rdname qcc-deprecated
#' @export
qcc.groups  <- function(...) 
{
  .Deprecated("qccGroups", package = "qcc")
  qccGroups(...)
}
