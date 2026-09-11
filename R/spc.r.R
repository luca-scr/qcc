#' Statistics used in computing and drawing a Shewhart R chart
#'
#' These functions are used to compute statistics required by the R chart.
#'
#' @inheritParams spc_common data center nsigmas conf
#' @param sizes samples sizes. Optional
#' @param std.dev within group standard deviation. Optional for `sd.R`
#'   function, required for `limits.R`. See [sd.xbar()].
#' @param ... catches further ignored arguments.
#' @return The function `stats.R` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.R` returns `std.dev` the standard deviation of
#' the statistic charted.
#'
#' The function `limits.R` returns a matrix with lower and upper control
#' limits.
#' @inherit spc_common author seealso
#' @references `r refs("montgomery2013", "wetherill_brown_1991")`
#' @name stats.R
NULL

#' @rdname stats.R
#' @export
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

#' @rdname stats.R
#' @export
sd.R <- function(data, sizes, std.dev = c("UWAVE-R", "MVLUE-R"), ...)
{
  if (!is.numeric(std.dev))
     std.dev <- match.arg(std.dev)
  sd.xbar(data, sizes, std.dev)
}

#' @rdname stats.R
#' @export
limits.R <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  if (length(unique(sizes))==1) sizes <- sizes[1]
  if (is.null(conf)) 
     {
       se.R <- .d3(sizes) * std.dev
       lcl <- pmax(0, center - nsigmas * se.R)
       ucl <- center + nsigmas * se.R
     }
  else 
     { if (conf > 0 && conf < 1) 
          # FIX: replace qtukey with a more precise implementation?
          { ucl <- qtukey(1 - (1 - conf)/2, sizes, 1e100) * std.dev
            lcl <- qtukey((1 - conf)/2, sizes, 1e100) * std.dev
          }
       else stop("invalid conf argument. See help.")
     }
  new_limits(lcl,ucl)
}
