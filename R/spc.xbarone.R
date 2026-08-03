#' Statistics used in computing and drawing a Shewhart xbar chart for one-at-time data
#'
#' These functions are used to compute statistics required by the xbar chart
#' for one-at-time data.
#'
#' Methods available for estimating the process standard deviation:
#'
#' - `"MR"`: moving range; this estimate is based on the scaled mean of moving
#'   ranges.
#' - `"SD"`: sample standard deviation; this estimate is defined as
#'   `sd(x) / cd(n)`, where `n` is the number of individual measurements of
#'   `x`.
#'
#' @inheritParams spc_common data center nsigmas conf
#' @param sizes samples sizes. Not needed, `size = 1` is used.
#' @param r number of successive pairs of observations for computing the
#' standard deviation based on moving ranges of r points.
#' @param std.dev within group standard deviation. Optional for
#'   `sd.xbar.one` function, required for `limits.xbar.one`. See details.
#' @param ... catches further ignored arguments.
#' @return The function `stats.xbar.one` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.xbar.one` returns `std.dev` the standard
#' deviation of the statistic charted.
#'
#' The function `limits.xbar.one` returns a matrix with lower and upper
#' control limits.
#' @inherit spc_common author seealso
#' @references `r refs("montgomery2013", "ryan_2011", "wetherill_brown_1991")`
#' @name stats.xbar.one
#' @examples
#' x <- antifreeze[["water"]] # See `?antifreeze`
#' # 1) using MR (default)
#' qcc(x, type="xbar.one", data.name="Water content (in ppm) of batches of antifreeze")
#' # 2) using SD
#' qcc(x, type="xbar.one", std.dev = "SD", data.name="Water content (in ppm) of batches of antifreeze")
#'
#' # "as the size increases further, we would expect sigma-hat to settle down
#' #  at a value close to the overall sigma-hat" (Wetherill and Brown, 1991, p. 121)
#' sigma <- NA
#' k <- 2:24
#' for (j in k) sigma[j] <- sd.xbar.one(x, k=j)
#' plot(k, sigma[k], type="b")     # plot estimates of sigma for 
#' abline(h=sd(x), col=2, lty=2)   # different values of k
NULL

# TODO: Estimator efficiency should be discussed in a vignette. Currently at docs @examples

#' @rdname stats.xbar.one
#' @export
stats.xbar.one <- function(data, sizes)
{
  statistics <- as.vector(data)
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

#' @rdname stats.xbar.one
#' @export
# PERF: Replace apply call with matrixStats
sd.xbar.one <- function(data, sizes, std.dev = c("MR", "SD"), r = 2, ...) {
  data <- as.vector(data)
  if(is.numeric(std.dev))
    return(std.dev)

  std.dev <- match.arg(std.dev)
  switch(std.dev, 
    "MR" = {
      windows <- embed(data[!is.na(data)], r)
      moving_ranges <- apply(windows, 1L, function(x) {
        diff(range(x))
      })
      mean(moving_ranges) / .d2(r)
    },
    "SD" = { 
      sd(data, na.rm = TRUE)/.c4(sum(!is.na(data)))
    })
}


#' @rdname stats.xbar.one
#' @export
limits.xbar.one <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  se.stats <- std.dev

  if (!is.null(conf)) {
    if (!is.numeric(conf) || length(conf) != 1L || conf <= 0 || conf >= 1) {
      stop("invalid 'conf' argument. See help.")
    }

    nsigmas <- qnorm(1 - (1 - conf) / 2)
  }

  delta <- nsigmas * se.stats
  lcl <- center - delta
  ucl <- center + delta
  new_limits(lcl,ucl)
}
