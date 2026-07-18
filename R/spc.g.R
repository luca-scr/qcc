#' Statistics used in computing and drawing a Shewhart g chart
#'
#' These functions are used to compute statistics required by the g chart
#' (geometric distribution) for use with the qcc package.
#'
#' The g chart plots the number of non-events between events.  np charts do not
#' work well when the probability of an event is rare (see example below).
#' Instead of plotting the number of events, the g chart plots the number of
#' non-events between events.
#'
#' @aliases stats.g sd.g limits.g
#' @export stats.g
#' @export sd.g
#' @export limits.g
#' @param data the observed data values
#' @param center sample center statistic
#' @param sizes sample sizes (not used)
#' @param std.dev standard deviation of geometric distribution
#' @param nsigmas a numeric value specifying the number of sigmas to use for
#' computing control limits. It is ignored when the `conf` argument is
#' provided.
#' @param conf a numeric value in \eqn{(0,1)} specifying the confidence level
#' to use for computing control limits.
#' @param ... catches further ignored arguments.
#' @return The function `stats.g()` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.g()` returns `std.dev` the standard deviation
#' \eqn{sqrt(1-p)/p}.
#'
#' The function `limits.g()` returns a matrix with lower and upper control
#' limits.
#' @note The geometric distribution is quite skewed so it is best to set conf
#' at the required confidence interval (0 < conf < 1) rather than as a
#' multiplier of sigma.
#' @author Greg Snow (greg.snow@ihc.com)
#' @seealso [qcc()]
#' @references Kaminsky, FC et. al. (1992) *Statistical Control Charts
#' Based on a Geometric Distribution*, Journal of Quality Technology, 24, pp
#' 63--69.
#'
#' Yang, Z et. al. (2002) On the Performance of Geometric Charts with Estimated
#' Control Limits, *Journal of Quality Technology*, 34, pp 448--458.
#' @keywords hplot
#' @examples
#'
#' success  = rbinom(1000, 1, 0.01)
#' num.noevent  = diff(which(c(1,success)==1))-1
#' qcc(success, type = "np", sizes = 1)
#' qcc(num.noevent, type = "g")
#'
stats.g <- function (data, sizes) 
{
  statistics <- as.vector(data)
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

sd.g <- function (data, sizes, ...)
{
  data <- as.vector(data)
  p <- 1/mean(data)
  std.dev <- sqrt( (1-p) )/p
  return(std.dev)
}

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
  .construct_limits(lcl,ucl)
}
