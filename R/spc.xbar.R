#' Statistics used in computing and drawing a Shewhart xbar chart
#'
#' These functions are used to compute statistics required by the xbar chart.
#'
#' The following methods are available for estimating the process standard
#' deviation:
#'
#' - `"UWAVE-R"`: UnWeighted AVErage of within-group estimates based on
#'   within-group Ranges.
#' - `"UWAVE-SD"`: UnWeighted AVErage of within-group estimates based on
#'   within-group Standard Deviations.
#' - `"MVLUE-R"`: Minimum Variance Linear Unbiased Estimator computed as a
#'   weighted average of within-group estimates based on within-group Ranges.
#' - `"MVLUE-SD"`: Minimum Variance Linear Unbiased Estimator computed as a
#'   weighted average of within-group estimates based on within-group Standard
#'   Deviations.
#' - `"RMSDF"`: Root-Mean-Square estimator computed as a weighted average of
#'   within-group estimates based on within-group Standard Deviations.
#'
#' Depending on the chart, a method may be available or not, or set as the
#' default according to the following table:
#'
#' | Method | `"xbar"` | `"R"` | `"S"` |
#' | --- | --- | --- | --- |
#' | `"UWAVE-R"` | default | default | not available |
#' | `"UWAVE-SD"` | available | not available | default |
#' | `"MVLUE-R"` | available | available | not available |
#' | `"MVLUE-SD"` | available | not available | available |
#' | `"RMSDF"` | available | not available | available |
#'
#' Detailed definitions of formulae implemented are available in the SAS/QC
#' User's Guide.
#'
#' @inheritParams spc_common data center nsigmas conf
#' @param sizes samples sizes. Optional
#' @param std.dev within group standard deviation. Optional for `sd.xbar`
#'   function, required for `limits.xbar`. See details.
#' @param ... catches further ignored arguments.
#' @return The function `stats.xbar` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.xbar` returns `std.dev` the standard deviation of
#' the statistic charted. This is based on results from Burr (1969).
#'
#' The function `limits.xbar` returns a matrix with lower and upper
#' control limits.
#' @inherit spc_common author seealso
#' @references `r refs("burr_1969", "montgomery2013", "wetherill_brown_1991")`
#' @name stats.xbar
NULL

#' @rdname stats.xbar
#' @export
stats.xbar <- function(data, sizes)
{
  data <- as.matrix(data)
  if(missing(sizes))
    sizes <- .rowNobs(data, useNames = FALSE)
  statistics <- rowMeans(data, na.rm = TRUE)
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

#' @rdname stats.xbar
#' @export
sd.xbar <- function(data, sizes, std.dev = c("UWAVE-R", "UWAVE-SD", "MVLUE-R", "MVLUE-SD", "RMSDF"), ...)
{
  data <- as.matrix(data)
  if(missing(sizes))
    sizes <- .rowNobs(data, useNames = FALSE)
  if(any(sizes == 1))
    stop("group sizes must be larger than one")
  if(is.numeric(std.dev))
    return(std.dev)

  std.dev <- match.arg(std.dev)
  switch(std.dev,
         "UWAVE-R" = {
           R <- .rowRanges(data, na.rm = TRUE)
           sum(R/.d2(sizes))/length(sizes)
         },
         "UWAVE-SD" = {
           S <- rowSds(data, na.rm = TRUE)
           sum(S/.c4(sizes))/length(sizes)
         },
         "MVLUE-R" = {
           R <- .rowRanges(data, na.rm = TRUE)
           d2 <- .d2(sizes)
           w <- (d2/.d3(sizes))^2
           sum(R/d2*w)/sum(w)
         },
         "MVLUE-SD" = {
           S <- rowSds(data, na.rm = TRUE)
           c4 <- .c4(sizes)
           w <- c4^2 / (1 - c4^2)
           sum(S / c4 * w) / sum(w)
         },
         "RMSDF" = {
           S <- rowSds(data, na.rm = TRUE)
           w <- sizes - 1
           sqrt(sum(S^2 * w) / sum(w)) / .c4(sum(w) + 1)
         })
}

#' @rdname stats.xbar
#' @export
limits.xbar <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")

  if (length(unique(sizes))==1) sizes <- sizes[1]
  se.stats <- std.dev/sqrt(sizes)

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
