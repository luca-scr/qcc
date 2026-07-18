#' Statistics used in computing and drawing a Shewhart u chart
#'
#' These functions are used to compute statistics required by the u chart.
#'
#' @inheritParams spc_common data center sizes std.dev nsigmas conf
#' @param ... catches further ignored arguments.
#' @return The function `stats.u` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.u` returns `std.dev` the standard deviation of
#' the statistic charted.
#'
#' The function `limits.u` returns a matrix with lower and upper control
#' limits.
#' @inherit spc_common author seealso
#' @references `r refs("montgomery2013", "wetherill_brown_1991")`
#' @name stats.u
NULL

#' @rdname stats.u
#' @export
stats.u <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  statistics <- data/sizes
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

#' @rdname stats.u
#' @export
sd.u <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  std.dev <- sqrt(sum(data)/sum(sizes))
  return(std.dev)
}

#' @rdname stats.u
#' @export
limits.u <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  sizes <- as.vector(sizes)
  if (length(unique(sizes))==1) sizes <- sizes[1]
  limits.c(center * sizes, std.dev, sizes, nsigmas, conf) / sizes
}
