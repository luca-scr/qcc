#' Statistics used in computing and drawing a Shewhart np chart
#'
#' These functions are used to compute statistics required by the np chart.
#'
#' @inheritParams spc_common data center sizes std.dev nsigmas conf
#' @param ... catches further ignored arguments.
#' @return The function `stats.np` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.np` returns `std.dev` the standard deviation of
#' the statistic charted.
#'
#' The function `limits.np` returns a matrix with lower and upper control
#' limits.
#' @inherit spc_common author seealso
#' @references `r refs("montgomery2013", "wetherill_brown_1991")`
#' @name stats.np
NULL

#' @rdname stats.np
#' @export
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

#' @rdname stats.np
#' @export
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

#' @rdname stats.np
#' @export
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
  new_limits(lcl,ucl)
}
