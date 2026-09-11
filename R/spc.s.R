#' Functions to plot Shewhart S chart
#'
#' These functions are used to compute statistics required by the S chart.
#'
#'
#' @inheritParams spc_common data center nsigmas conf
#' @param sizes samples sizes. Optional
#' @param std.dev within group standard deviation. Optional for `sd.S`
#'   function, required for `limits.S`. See [sd.xbar()].
#' @param ... catches further ignored arguments.
#' @return The function `stats.S` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.S` returns `std.dev` the standard deviation of
#' the statistic charted.
#'
#' The function `limits.S` returns a matrix with lower and upper control
#' limits.
#' @inherit spc_common author seealso
#' @references `r refs("montgomery2013", "wetherill_brown_1991")`
#' @name stats.S
NULL

#' @rdname stats.S
#' @export
stats.S <- function(data, sizes)
{
  data <- as.matrix(data)
  if (missing(sizes))
     sizes <- .rowNobs(data, useNames = FALSE)
  if(ncol(data)==1) 
    { statistics <- as.vector(data) }
  else 
    { statistics <- .rowSds(data, na.rm = TRUE) }
  if (length(sizes == 1))
     sizes <- rep(sizes, length(statistics))
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

#' @rdname stats.S
#' @export
sd.S <- function(data, sizes, std.dev = c("UWAVE-SD", "MVLUE-SD", "RMSDF"), ...)
{
  if (!is.numeric(std.dev))
     std.dev <- match.arg(std.dev)
  sd.xbar(data, sizes, std.dev)
}

#' @rdname stats.S
#' @export
limits.S <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  if(length(unique(sizes))==1) sizes <- sizes[1]
  se.stats <- std.dev * sqrt(1 - .c4(sizes)^2)
  if (is.null(conf)) 
     { lcl <- pmax(0, center - nsigmas * se.stats)
       ucl <- center + nsigmas * se.stats
     }
  else 
     { if (conf > 0 & conf < 1) 
          { ucl <- std.dev * sqrt(qchisq(1 - (1 - conf)/2, sizes - 1)/
                                  (sizes - 1))
            lcl <- std.dev * sqrt(qchisq((1 - conf)/2, sizes - 1)/
                                  (sizes - 1))
          }
          else stop("invalid conf argument. See help.")
     }
  new_limits(lcl,ucl)
}
