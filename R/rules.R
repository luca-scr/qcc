# TODO: Refactor WER2, NEL7, NEL8

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
  if(!(inherits(object, "qcc") | inherits(object, "mqcc")))
    stop("input object must be of class 'qcc' or 'mqcc'")
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

qccRulesViolatingWER1 <- function(object, limits = object$limits)
{
  # Return cases beyond control limits (WER #1)
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

qccRulesViolatingWER4 <- function(object) qccRulesViolatingNEL2(object, run.length = 8)

# Nelson rules
#
# A process is out of control if any of the following occur:
# 1. One point plots outside 3-sigma control limits.
# 2. Nine points in a row plot on the same side of the center line.
# 3. Six points in a row are steadily increasing or decreasing.
# 4. Fourteen points in a row alternate up and down.
# 5. Two of three consecutive points plot beyond a 2-sigma limit.
# 6. Four of five consecutive points plot beyond a 1-sigma limit.
# 7. Fifteen points in a row plot within 1 sigma of the center line.
# 8. Eight points in a row plot outside 1 sigma on both sides of the center line.

qccRulesViolatingNEL1 <- function(object) qccRulesViolatingWER1(object, object$limits)

qccRulesViolatingNEL2 <- function(object, run.length = 9)
{
  # Return indices of points violating nine-point runs (Nelson #2)
  center <- object$center
  statistics <- c(object$statistics, object$newstats)
  diffs <- statistics - center
  viol.above <- qccRulesViolatingRun(diffs > 0, run.length)
  viol.below <- qccRulesViolatingRun(diffs < 0, run.length)
  return(c(viol.above, viol.below))
}

qccRulesViolatingNEL3 <- function(object)
{
  # Return indices of points violating trend runs (Nelson #3)
  run.length <- 6
  statistics <- c(object$statistics, object$newstats)
  diffs <- diff(statistics)
  viol.increase <- qccRulesViolatingRun(diffs > 0, run.length - 1) + 1
  viol.decrease <- qccRulesViolatingRun(diffs < 0, run.length - 1) + 1
  return(c(viol.increase, viol.decrease))
}

qccRulesViolatingNEL4 <- function(object)
{
  # Return indices of points violating alternating runs (Nelson #4)
  run.length <- 14
  statistics <- c(object$statistics, object$newstats)
  diffs <- diff(statistics)
  if(length(diffs) < 2) # Need at least 2 to detect alternating direction
    return(numeric())
  alternating <- diffs[-1] * diffs[-length(diffs)] < 0 # Identify where consecutive differences change sign
  violators <- qccRulesViolatingRun(alternating, run.length - 2) + 2
  return(violators)
}

qccRulesViolatingNEL5 <- function(object) qccRulesViolatingWER2(object)
qccRulesViolatingNEL6 <- function(object) qccRulesViolatingWER3(object)

qccRulesViolatingNEL7 <- function(object)
{
  # Return indices of points inside one-sigma limits (Nelson #7)
  run.length <- 15
  statistics <- c(object$statistics, object$newstats)
  limits <- qccRulesSigmaLimits(object, k = 1)
  inside <- statistics >= limits[,1] & statistics <= limits[,2]
  return(qccRulesViolatingRun(inside, run.length))
}

qccRulesViolatingNEL8 <- function(object)
{
  # Return indices of points outside one-sigma limits on both sides (Nelson #8)
  run.length <- 8
  statistics <- c(object$statistics, object$newstats)
  limits <- qccRulesSigmaLimits(object, k = 1)
  if(length(statistics) < run.length)
    return(numeric())
  above <- statistics > limits[,2]
  below <- statistics < limits[,1]
  outside <- above | below
  outside.windows <- embed(outside, run.length)
  above.windows <- embed(above, run.length)
  below.windows <- embed(below, run.length)
  violators <- which(apply(outside.windows, 1, all) &
                     apply(above.windows, 1, any) &
                     apply(below.windows, 1, any))
  return(violators + run.length - 1)
}

qccRulesViolatingRun <- function(condition, run.length) {
  # Returns the indices of elements that are part of a long enough run of TRUE values
  # only starts counting from the point where the run first reaches the required length.
  if(!length(condition))
    return(numeric())

  condition[is.na(condition)] <- FALSE
  runs <- rle(condition)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  violating.runs <- which(runs$values & runs$lengths >= run.length)
  violators <- numeric()
  if (length(violating.runs)) {
    for(i in violating.runs)
    violators <- c(violators, (starts[i] + run.length - 1):ends[i])
  }
  return(violators)
}

qccRulesSigmaLimits <- function(object, k)
{
  statistics <- c(object$statistics, object$newstats)
  limits <- paste("limits.", object$type, sep = "")
  limits <- do.call(limits, list(center = object$center,
                                 std.dev = object$std.dev,
                                 sizes = c(object$sizes, object$newsizes),
                                 nsigmas = object$nsigmas*k/3))
  if(nrow(limits) == 1)
    limits <- limits[rep(1, length(statistics)),, drop = FALSE]
  return(limits)
}
