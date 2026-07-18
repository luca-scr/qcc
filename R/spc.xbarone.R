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
#' @aliases stats.xbar.one sd.xbar.one limits.xbar.one
#' @export stats.xbar.one
#' @export sd.xbar.one
#' @export limits.xbar.one
#' @param data the observed data values
#' @param center sample/group center statistic.
#' @param sizes samples sizes. Not needed, `size = 1` is used.
#' @param r number of successive pairs of observations for computing the
#' standard deviation based on moving ranges of r points.
#' @param std.dev within group standard deviation. Optional for
#' `sd.xbar.one` function, required for `limits.xbar.one`. See
#' details.
#' @param nsigmas a numeric value specifying the number of sigmas to use for
#' computing control limits. It is ignored when the `conf` argument is
#' provided.
#' @param conf a numeric value in \eqn{(0,1)} specifying the confidence level
#' to use for computing control limits.
#' @param ... catches further ignored arguments.
#' @return The function `stats.xbar.one` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.xbar.one` returns `std.dev` the standard
#' deviation of the statistic charted.
#'
#' The function `limits.xbar.one` returns a matrix with lower and upper
#' control limits.
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @references Montgomery, D.C. (2013) *Introduction to Statistical
#' Quality Control*, 7th ed. New York: John Wiley & Sons.
#'
#' Ryan, T. P. (2011), *Statistical Methods for Quality Improvement*, 3rd
#' ed. New York: John Wiley & Sons, Inc.
#'
#' Wetherill, G.B. and Brown, D.W. (1991) *Statistical Process Control*.
#' New York: Chapman & Hall.
#' @keywords htest hplot
#' @examples
#'
#' # Water content of antifreeze data (Wetherill and Brown, 1991, p. 120)
#' x  = c(2.23, 2.53, 2.62, 2.63, 2.58, 2.44, 2.49, 2.34, 2.95, 2.54, 2.60, 2.45,
#'        2.17, 2.58, 2.57, 2.44, 2.38, 2.23, 2.23, 2.54, 2.66, 2.84, 2.81, 2.39,
#'        2.56, 2.70, 3.00, 2.81, 2.77, 2.89, 2.54, 2.98, 2.35, 2.53)
#' # the Shewhart control chart for one-at-time data
#' # 1) using MR (default)
#' qcc(x, type="xbar.one", data.name="Water content (in ppm) of batches of antifreeze")
#' # 2) using SD
#' qcc(x, type="xbar.one", std.dev = "SD", data.name="Water content (in ppm) of batches of antifreeze")
#'
#' # "as the size increases further, we would expect sigma-hat to settle down
#' #  at a value close to the overall sigma-hat" (Wetherill and Brown, 1991,
#' # p. 121)
#' sigma  = NA
#' k  = 2:24
#' for (j in k)
#'     sigma[j]  = sd.xbar.one(x, k=j)
#' plot(k, sigma[k], type="b")     # plot estimates of sigma for 
#' abline(h=sd(x), col=2, lty=2)   # different values of k
#'
stats.xbar.one <- function(data, sizes)
{
  statistics <- as.vector(data)
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

sd.xbar.one <- function(data, sizes, std.dev = c("MR", "SD"), r = 2, ...)
{
  data <- as.vector(data)
  n <- length(data)
  if(!is.numeric(std.dev)) 
     std.dev <- match.arg(std.dev)
  if(is.numeric(std.dev)) 
    { sd <- std.dev }
  else
    { switch(std.dev, 
             "MR" = {
                data <- data[!is.na(data)]
                d2 <- qcc.options("exp.R.unscaled")
                moving_ranges <- apply(embed(data, r), 1L, function(x) {
                  diff(range(x))
                })
                sd <- mean(moving_ranges) / d2[r]
             },
             "SD" = { sd <- sd(data, na.rm = TRUE)/qcc.c4(sum(!is.na(data))) },
             sd <- NULL)
    }
  return(sd)
}


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
  .construct_limits(lcl,ucl)
}
