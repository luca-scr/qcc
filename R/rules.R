# Western Electric rules 
#
# A process is out of control if either
# 1. One point plots outside 3-sigma control limits.
# 2. Two of three consecutive points plot beyond a 2-sigma limit.
# 3. Four of five consecutive points plot beyond a 1-sigma limit.
# 4. Eight consecutive points plot on one side of the center line.

qccRules <- function(object, rules = object$rules)
{
  # Return a vector of indices for cases (statistics & new.statistics) 
  # in object violating specified rules (NA if no rule is violated)
  rules <- as.numeric(rules)
  if(!inherits(object, "qcc"))
    stop("input object must be of class 'qcc'")
  stats <- c(object$statistics, object$newstats)
  out <- rep(NA, length(stats))
  if(any(rules == 4)) 
  {  
    wer <- qccRulesViolatingWER4(object)
    out[wer] <- 4
  }
  if(any(rules == 3)) 
  {  
    wer <- qccRulesViolatingWER3(object)
    out[wer] <- 3
  }
  if(any(rules == 2)) 
  {  
    wer <- qccRulesViolatingWER2(object)
    out[wer] <- 2
  }
  if(any(rules == 1)) 
  {  
    wer <- qccRulesViolatingWER1(object)
    out[wer] <- 1
  }
  attr(out, "WesternElectricRules") <- rules
  return(out)
}

qccRulesViolatingWER1 <- function(object)
{
  # Return cases beyond control limits (WER #1)
  limits <- object$limits
  statistics <- c(object$statistics, object$newstats) 
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.above.ucl <- seq(along = statistics)[statistics > ucl]
  index.below.lcl <- seq(along = statistics)[statistics < lcl]
  return(c(index.above.ucl, index.below.lcl))
}

qccRulesViolatingWER2 <- function(object, 
                                  run.points = 2,
                                  run.length = 3,
                                  k = object$nsigmas*2/3)
{
  # Return indices of points violating runs
  center     <- object$center
  statistics <- c(object$statistics, object$newstats)
  limits     <- paste("limits.", object$type, sep = "")
  limits     <- do.call(limits, list(center = object$center, 
                                     std.dev = object$std.dev,
                                     sizes = c(object$sizes, object$newsizes),
                                     nsigmas = k))
  i <- if(nrow(limits) > 1) seq(run.length, length(statistics)) else 1
  viol.above <- embed(statistics, run.length) > limits[i,2]
  viol.above <- which(apply(viol.above, 1, sum) >= run.points & viol.above[,1])
  viol.above <- viol.above + (run.length-1)
  viol.below <- embed(statistics, run.length) < limits[i,1]
  viol.below <- which(apply(viol.below, 1, sum) >= run.points & viol.below[,1])
  viol.below <- viol.below + (run.length-1)
  return(c(viol.above, viol.below))
}

qccRulesViolatingWER3 <- function(object, ...)
{
  qccRulesViolatingWER2(object, 
                        run.points = 4,
                        run.length = 5,
                        k = object$nsigmas*1/3)
}  

qccRulesViolatingWER4 <- function(object)
{
  # Return indices of points violating runs (WER #4)
  run.length <- 8
  center <- object$center
  statistics <- c(object$statistics, object$newstats)
  cl <- object$limits
  diffs <- statistics - center
  diffs[diffs > 0] <- 1
  diffs[diffs < 0] <- -1
  runs <- rle(diffs)
  vruns <- rep(runs$lengths >= run.length, runs$lengths)
  vruns.above <- (vruns & (diffs > 0))
  vruns.below <- (vruns & (diffs < 0))
  rvruns.above <- rle(vruns.above)
  rvruns.below <- rle(vruns.below)
  vbeg.above <- cumsum(rvruns.above$lengths)[rvruns.above$values] -
    (rvruns.above$lengths - run.length)[rvruns.above$values]
  vend.above <- cumsum(rvruns.above$lengths)[rvruns.above$values]
  vbeg.below <- cumsum(rvruns.below$lengths)[rvruns.below$values] -
    (rvruns.below$lengths - run.length)[rvruns.below$values]
  vend.below <- cumsum(rvruns.below$lengths)[rvruns.below$values]
  violators <- numeric()
  if (length(vbeg.above)) 
  { for (i in 1:length(vbeg.above))
    violators <- c(violators, vbeg.above[i]:vend.above[i]) }
  if (length(vbeg.below)) 
  { for (i in 1:length(vbeg.below))
    violators <- c(violators, vbeg.below[i]:vend.below[i]) }
  return(violators)
}



#
# Functions used to signal points out of control in 'qcc' version <= 2.x
#
# TODO: to remove when new rules are ok

shewhart.rules <- function(object, limits = object$limits, run.length = qcc.options("run.length"))
{
# Return a list of cases beyond limits and violating runs
  bl <- beyond.limits(object, limits = limits)
  vr <- violating.runs(object, run.length = run.length)
  list(beyond.limits = bl, violating.runs = vr)
}

# shewhartRules <- function(object, 
#                            limits = object$limits, 
#                            run.length = qcc.options("run.length"))
# {
#   # Return a vector of indices for cases (statistics & new.statistics) 
#   # violating the two basic Shewhart rules (NA if no rule is violated)
#   stats <- c(object$statistics, object$newstats)
#   out <- rep(NA, length(stats))
#   bl <- beyond.limits(object, limits = limits)
#   vr <- violating.runs(object, run.length = run.length)
#   out[vr] <- 2
#   out[bl] <- 1
#   attr(out, "rules") <- "shr"
#   return(out)
# }

beyond.limits <- function(object, limits = object$limits)
{
  # Return cases beyond limits
  statistics <- c(object$statistics, object$newstats) 
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.above.ucl <- seq(along = statistics)[statistics > ucl]
  index.below.lcl <- seq(along = statistics)[statistics < lcl]
  return(c(index.above.ucl, index.below.lcl))
}

violating.runs <- function(object, run.length = qcc.options("run.length"))
{
  # Return indices of points violating runs
  if(run.length == 0)
    return(numeric())
  center <- object$center
  statistics <- c(object$statistics, object$newstats)
  cl <- object$limits
  diffs <- statistics - center
  diffs[diffs > 0] <- 1
  diffs[diffs < 0] <- -1
  runs <- rle(diffs)
  vruns <- rep(runs$lengths >= run.length, runs$lengths)
  vruns.above <- (vruns & (diffs > 0))
  vruns.below <- (vruns & (diffs < 0))
  rvruns.above <- rle(vruns.above)
  rvruns.below <- rle(vruns.below)
  vbeg.above <- cumsum(rvruns.above$lengths)[rvruns.above$values] -
    (rvruns.above$lengths - run.length)[rvruns.above$values]
  vend.above <- cumsum(rvruns.above$lengths)[rvruns.above$values]
  vbeg.below <- cumsum(rvruns.below$lengths)[rvruns.below$values] -
    (rvruns.below$lengths - run.length)[rvruns.below$values]
  vend.below <- cumsum(rvruns.below$lengths)[rvruns.below$values]
  violators <- numeric()
  if (length(vbeg.above)) 
  { for (i in 1:length(vbeg.above))
    violators <- c(violators, vbeg.above[i]:vend.above[i]) }
  if (length(vbeg.below)) 
  { for (i in 1:length(vbeg.below))
    violators <- c(violators, vbeg.below[i]:vend.below[i]) }
  return(violators)
}

