"stats.g" <- function (data, sizes) 
{
  statistics <- as.vector(data)
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

"sd.g" <- function (data, sizes, ...)
{
  data <- as.vector(data)
  p <- 1/mean(data)
  std.dev <- sqrt( (1-p) )/p
  return(std.dev)
}

"limits.g" <- function (center, std.dev, sizes, conf) 
{
  if (conf >= 1) {
    p <- 1/center
    lcl <- center - conf * sqrt(1-p)/p
    lcl[lcl < 0] <- 0
    ucl <- center + conf * sqrt(1-p)/p
    warning("The Geometric distribution is quite skewed, it is better to set conf at the required confidence level (0 < conf < 1) instead of as a multiplier of sigma.")
  }
  else {
    if (conf > 0 & conf < 1) {
      p <- 1/center
      ucl <- qgeom(1 - (1 - conf)/2, p)
      lcl <- qgeom((1 - conf)/2, p)
    }
    else stop("invalid conf argument. See help.")
  }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

