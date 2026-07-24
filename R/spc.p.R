#' Statistics used in computing and drawing a Shewhart p chart
#'
#' These functions are used to compute statistics required by the p chart.
#'
#' @inheritParams spc_common data center sizes std.dev nsigmas conf
#' @param ... catches further ignored arguments.
#' @return The function `stats.p` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.p` returns `std.dev` the standard deviation of
#' the statistic charted.
#'
#' The function `limits.p` returns a matrix with lower and upper control
#' limits.
#' @inherit spc_common author seealso
#' @references `r refs("montgomery2013", "wetherill_brown_1991")`
#' @name stats.p
NULL

#' @rdname stats.p
#' @export
stats.p <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  list(statistics = data/sizes, center = pbar)
}

#' @rdname stats.p
#' @export
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

#' @rdname stats.p
#' @export
limits.p <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{ 
  limits.np(center * sizes, std.dev, sizes, nsigmas, conf) / sizes
}
