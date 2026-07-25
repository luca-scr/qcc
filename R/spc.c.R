#' Functions to plot Shewhart c chart
#'
#' Statistics used in computing and drawing a Shewhart c chart.
#'
#' @inheritParams spc_common data center sizes std.dev nsigmas conf
#' @param ... catches further ignored arguments.
#' @return The function `stats.c` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.c` returns `std.dev` the standard deviation of
#' the statistic charted.
#'
#' The function `limits.c` returns a matrix with lower and upper control
#' limits.
#' @inherit spc_common author seealso
#' @references `r refs("montgomery2013", "wetherill_brown_1991")`
#' @name stats.c
NULL

#' @rdname stats.c
#' @export
stats.c <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) != 1)
     stop("all sizes must be equal for a c chart")
  statistics <- data
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

#' @rdname stats.c
#' @export
sd.c <- function(data, sizes, ...)
{
  data <- as.vector(data)
  std.dev <- sqrt(mean(data))
  return(std.dev)
}

#' @rdname stats.c
#' @export
limits.c <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  if (is.null(conf))
     { lcl <- pmax(0, center - nsigmas * sqrt(center))
       ucl <- center + nsigmas * sqrt(center)
     }
  else 
     { if (conf > 0 & conf < 1) 
          { ucl <- qpois(1 - (1 - conf)/2, center)
            lcl <- qpois((1 - conf)/2, center)
          }
       else stop("invalid conf argument. See help.")
     }
  new_limits(lcl,ucl)
}
