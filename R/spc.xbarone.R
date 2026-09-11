#' Statistics used in computing and drawing a Shewhart xbar chart for one-at-time data
#'
#' These functions are used to compute statistics required by the xbar chart
#' for one-at-time data.
#'
#' Methods available for estimating the process standard deviation:
#'
#' - `"MR"`: moving range; this estimator is based on the scaled mean of moving
#'   ranges.
#' - `"MMR"`: median moving range; this estimator is defined as
#'   `median(moving_ranges) / d4(r)`, where [d4()] is the median range of `r`
#'   independent standard normal observations.
#' - `"SD"`: sample standard deviation; this estimator is defined as
#'   `sd(x) / c4(n)`, where `n` is the number of individual measurements of
#'   `x`.
#' - `"MSSD"`: mean squared successive differences; this estimator is defined as
#'   \eqn{\sqrt{\frac{1}{2(n-1)}\sum_{i=1}^{n-1}(x_{i+1}-x_i)^2}/c_4'(n)}{sqrt(mean(diff(x)^2) / 2) / c4'(n)}.
#'   The correction factor \eqn{c_4'(n)}{c4'(n)} differs from [c4()] by accounting
#'   for the dependence between successive differences.
#'
#' For independent normal observations from an in-control process, `"MSSD"`
#' generally has lower sampling variance than `"MR"` with \code{r = 2}. It can
#' also perform well under some changes in the process mean. However,
#' squaring the differences makes it sensitive to outliers. Serial correlation
#' in the observations can also bias the estimate; the correction factor
#' accounts for dependence between differences, not between observations.
#' See Braun and Park (2008) for comparisons of estimators for individuals charts.
#'
#' The current implementation of the MMR estimator is not corrected for
#' finite-sample bias.
#'
#' Missing values are omitted when estimating the standard deviation. For
#' `"MR"`, `"MMR"`, and `"MSSD"`, successive observations refer to the remaining values
#' in their original order, including pairs across gaps left by missing values.
#'
#' @inheritParams spc_common data center nsigmas conf
#' @param sizes samples sizes. Not needed, `size = 1` is used.
#' @param r number of successive observations in each moving range for the
#'   `"MR"` and `"MMR"` estimators. Ignored for `"SD"` and `"MSSD"`.
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
#' @references `r refs("montgomery2013", "ryan_2011", "wetherill_brown_1991", "braun_park_2008")`
#' @name stats.xbar.one
#' @examples
#' x <- antifreeze[["water"]] # See `?antifreeze`
#' # 1) using MR (default)
#' qcc(x, type="xbar.one", data.name="Water content (in ppm) of batches of antifreeze")
#' # 2) using SD
#' qcc(x, type="xbar.one", std.dev = "SD", data.name="Water content (in ppm) of batches of antifreeze")
#' # 3) using bias-corrected MSSD
#' qcc(x, type="xbar.one", std.dev = "MSSD", data.name="Water content (in ppm) of batches of antifreeze")
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
# TODO: 2. Implement the median absolute deviation estimator (MAD)
# TODO: 3. implement the Ciminera-Tukey Estimator
# TODO: 4. implement the Boyles Dynamic Linear model estimator (DLM)
# TODO: explain in docs why MMR is bad when should you use it.
sd.xbar.one <- function(data, sizes, std.dev = c("MR", "SD", "MSSD", "MMR"), r = 2, ...) {
  data <- as.vector(data)
  if(is.numeric(std.dev))
    return(std.dev)

  std.dev <- match.arg(std.dev)
  switch(std.dev, 
    "MR" =,
    "MMR" = {
      windows <- embed(data[!is.na(data)], r)
      moving_ranges <- .rowRanges(windows)
      if (std.dev == "MMR")
        median(moving_ranges) / .d4(r)
      else
        mean(moving_ranges) / .d2(r)
    },
    "SD" = { 
      sd(data, na.rm = TRUE)/.c4(sum(!is.na(data)))
    },
    "MSSD" = {
      data <- data[!is.na(data)]
      n <- length(data)
      if (n < 2L)
        return(NA_real_)

      sqrt(mean(diff(data)^2) / 2) / .c4_mssd(n)
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
