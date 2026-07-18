#' Statistics used in computing and drawing a Shewhart R chart
#'
#' These functions are used to compute statistics required by the R chart.
#'
#'
#' @aliases stats.R sd.R limits.R
#' @export stats.R
#' @export sd.R
#' @export limits.R
#' @param data the observed data values
#' @param center sample/group center statistic.
#' @param sizes samples sizes. Optional
#' @param std.dev within group standard deviation. Optional for `sd.R`
#' function, required for `limits.R`. See [sd.xbar()].
#' @param nsigmas a numeric value specifying the number of sigmas to use for
#' computing control limits. It is ignored when the `conf` argument is
#' provided.
#' @param conf a numeric value in \eqn{(0,1)} specifying the confidence level
#' to use for computing control limits.
#' @param ... catches further ignored arguments.
#' @return The function `stats.R` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.R` returns `std.dev` the standard deviation of
#' the statistic charted.
#'
#' The function `limits.R` returns a matrix with lower and upper control
#' limits.
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @references Montgomery, D.C. (2013) *Introduction to Statistical
#' Quality Control*, 7th ed. New York: John Wiley & Sons.
#'
#' Wetherill, G.B. and Brown, D.W. (1991) *Statistical Process Control*.
#' New York: Chapman & Hall.
#' @keywords htest hplot
stats.R <- function(data, sizes)
{
  data <- as.matrix(data)
  if (missing(sizes))
     sizes <- as.integer(rowSums(!is.na(data)))
  if(ncol(data)==1) 
    { statistics <- as.vector(data) }
  else 
    { statistics <- apply(data, 1, function(x) diff(range(x, na.rm=TRUE))) }
  if (length(sizes == 1))
     sizes <- rep(sizes, length(statistics))
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.R <- function(data, sizes, std.dev = c("UWAVE-R", "MVLUE-R"), ...)
{
  if (!is.numeric(std.dev))
     std.dev <- match.arg(std.dev)
  sd.xbar(data, sizes, std.dev)
}

limits.R <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  if (length(unique(sizes))==1) sizes <- sizes[1]
  se.R.unscaled <- qcc.options("se.R.unscaled")
  Rtab <- length(se.R.unscaled)
  if (is.null(conf)) 
     { if (any(sizes > Rtab))
          stop(paste("group size must be less than", 
                      Rtab + 1, "when giving nsigmas"))
       se.R <- se.R.unscaled[sizes] * std.dev
       lcl <- pmax(0, center - nsigmas * se.R)
       ucl <- center + nsigmas * se.R
     }
  else 
     { if (conf > 0 && conf < 1) 
          { ucl <- qtukey(1 - (1 - conf)/2, sizes, 1e100) * std.dev
            lcl <- qtukey((1 - conf)/2, sizes, 1e100) * std.dev
          }
       else stop("invalid conf argument. See help.")
     }
  .construct_limits(lcl,ucl)
}
