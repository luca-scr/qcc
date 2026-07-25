#' Statistics used in computing and drawing a Shewhart g chart
#'
#' These functions are used to compute statistics required by the g chart
#' (geometric distribution) for use with the qcc package.
#'
#' The g chart plots the number of non-events between events. np charts do not
#' work well when the probability of an event is rare (see example below).
#' Instead of plotting the number of events, the g chart plots the number of
#' non-events between events.
#'
#' The geometric distribution is quite skewed so it is best to set `conf`
#' at the required confidence interval (0 < conf < 1) rather than as a
#' multiplier of sigma.
#'
#' @inheritParams spc_common data center nsigmas conf
#' @param sizes sample sizes (not used)
#' @param std.dev standard deviation of geometric distribution
#' @param ... catches further ignored arguments.
#' @return The function `stats.g()` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.g()` returns `std.dev` the standard deviation
#' \eqn{sqrt(1-p)/p}.
#'
#' The function `limits.g()` returns a matrix with lower and upper control
#' limits.
#' @author Greg Snow (greg.snow@ihc.com)
#' @inherit spc_common seealso
#' @references `r refs("kaminsky_1992", "yang_2002")`
#' @name stats.g
#' @examples
#'
#' success  = rbinom(1000, 1, 0.01)
#' num.noevent  = diff(which(c(1,success)==1))-1
#' qcc(success, type = "np", sizes = 1)
#' qcc(num.noevent, type = "g")
NULL

#' @rdname stats.g
#' @export
stats.g <- function (data, sizes) 
{
  statistics <- as.vector(data)
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

#' @rdname stats.g
#' @export
sd.g <- function (data, sizes, ...)
{
  data <- as.vector(data)
  p <- 1/mean(data)
  std.dev <- sqrt( (1-p) )/p
  return(std.dev)
}

#' @rdname stats.g
#' @export
limits.g <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  if (is.null(conf)) 
    {
      p <- 1/center
      lcl <- pmax(0, center - conf * sqrt(1-p)/p)
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
  new_limits(lcl,ucl)
}
