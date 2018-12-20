#
#  Deprecated functions
#

pareto.chart <- function(...) 
{
  .Deprecated("paretoChart", package = "qcc")
  paretoChart(...)
}

process.capability <- function(...) 
{
  .Deprecated("processCapability", package = "qcc")
  processCapability(...)
}

oc.curves <- function(...) 
{
  .Deprecated("ocCurves", package = "qcc")
  ocCurves(...)
}

qcc.overdispersion.test <- function(...) 
{
  .Deprecated("qccOverdispersionTest", package = "qcc")
  qccOverdispersionTest(...)
}

cause.and.effect <- function(...) 
{
  .Deprecated("causeEffectDiagram", package = "qcc")
  causeEffectDiagram(...)
}

qcc.groups  <- function(...) 
{
  .Deprecated("qccGroups", package = "qcc")
  qccGroups(...)
}
