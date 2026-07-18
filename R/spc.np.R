#' Statistics used in computing and drawing a Shewhart np chart
#'
#' These functions are used to compute statistics required by the np chart.
#'
#'
#' @aliases stats.np sd.np limits.np
#' @export stats.np
#' @export sd.np
#' @export limits.np
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
#' @return The function `stats.np` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.np` returns `std.dev` the standard deviation of
#' the statistic charted.
#'
#' The function `limits.np` returns a matrix with lower and upper control
#' limits.
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @references Montgomery, D.C. (2013) *Introduction to Statistical
#' Quality Control*, 7th ed. New York: John Wiley & Sons.
#'
#' Wetherill, G.B. and Brown, D.W. (1991) *Statistical Process Control*.
#' New York: Chapman & Hall.
#' @keywords htest hplot
stats.np <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  center <- sizes * pbar
  if (length(unique(center)) == 1)
     center <- center[1]
  list(statistics = data, center = center)
}

sd.np <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  std.dev <- sqrt(sizes * pbar * (1 - pbar))
  if (length(unique(std.dev)) == 1)
     std.dev <- std.dev[1]
  return(std.dev)
}

limits.np <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{ 
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) == 1) sizes <- sizes[1]
  pbar <- mean(center / sizes)
  if (is.null(conf))
     { tol <- nsigmas * sqrt(pbar * (1 - pbar) * sizes)
       lcl <- pmax(center - tol, 0)
       ucl <- pmin(center + tol, sizes)
     }
  else
     { if (conf > 0 & conf < 1)
          { lcl <- qbinom((1 - conf)/2, sizes, pbar)
            ucl <- qbinom((1 - conf)/2, sizes, pbar, lower.tail = FALSE)
          }
       else stop("invalid conf argument. See help.")
     }
  .construct_limits(lcl,ucl)
}
