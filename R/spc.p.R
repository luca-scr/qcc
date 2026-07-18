#' Statistics used in computing and drawing a Shewhart p chart
#'
#' These functions are used to compute statistics required by the p chart.
#'
#'
#' @aliases stats.p sd.p limits.p
#' @export stats.p
#' @export sd.p
#' @export limits.p
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
#' @return The function `stats.p` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.p` returns `std.dev` the standard deviation of
#' the statistic charted.
#'
#' The function `limits.p` returns a matrix with lower and upper control
#' limits.
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @references Montgomery, D.C. (2013) *Introduction to Statistical
#' Quality Control*, 7th ed. New York: John Wiley & Sons.
#'
#' Wetherill, G.B. and Brown, D.W. (1991) *Statistical Process Control*.
#' New York: Chapman & Hall.
#' @keywords htest hplot
stats.p <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  list(statistics = data/sizes, center = pbar)
}

sd.p <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  std.dev <- sqrt(pbar * (1 - pbar) / sizes)
  if (length(unique(std.dev)) == 1)
     std.dev <- std.dev[1]
  return(std.dev)
}

limits.p <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{ 
  limits.np(center * sizes, std.dev, sizes, nsigmas, conf) / sizes
}
