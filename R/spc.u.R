#' Statistics used in computing and drawing a Shewhart u chart
#'
#' These functions are used to compute statistics required by the u chart.
#'
#'
#' @aliases stats.u sd.u limits.u
#' @export stats.u
#' @export sd.u
#' @export limits.u
#' @param data the observed data values
#' @param center sample/group center statistic.
#' @param sizes samples sizes.
#' @param std.dev within group standard deviation.
#' @param nsigmas a numeric value specifying the number of sigmas to use for
#' computing control limits. It is ignored when the `conf` argument is
#' provided.
#' @param conf a numeric value in \eqn{(0,1)} specifying the confidence level
#' to use for computing control limits.
#' @param ... catches further ignored arguments.
#' @return The function `stats.u` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.u` returns `std.dev` the standard deviation of
#' the statistic charted.
#'
#' The function `limits.u` returns a matrix with lower and upper control
#' limits.
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @references Montgomery, D.C. (2013) *Introduction to Statistical
#' Quality Control*, 7th ed. New York: John Wiley & Sons.
#'
#' Wetherill, G.B. and Brown, D.W. (1991) *Statistical Process Control*.
#' New York: Chapman & Hall.
#' @keywords htest hplot
stats.u <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  statistics <- data/sizes
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.u <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  std.dev <- sqrt(sum(data)/sum(sizes))
  return(std.dev)
}

limits.u <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  sizes <- as.vector(sizes)
  if (length(unique(sizes))==1) sizes <- sizes[1]
  limits.c(center * sizes, std.dev, sizes, nsigmas, conf) / sizes
}
