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
#' @aliases stats.xbar sd.xbar limits.xbar
#' @export stats.xbar
#' @export sd.xbar
#' @export limits.xbar
#' @param data the observed data values
#' @param center sample/group center statistic
#' @param sizes samples sizes. Optional
#' @param std.dev within group standard deviation. Optional for `sd.xbar`
#' function, required for `limits.xbar`. See details.
#' @param nsigmas a numeric value specifying the number of sigmas to use for
#' computing control limits. It is ignored when the `conf` argument is
#' provided.
#' @param conf a numeric value in \eqn{(0,1)} specifying the confidence level
#' to use for computing control limits.
#' @param ... catches further ignored arguments.
#' @return The function `stats.xbar` returns a list with components
#' `statistics` and `center`.
#'
#' The function `sd.xbar` returns `std.dev` the standard deviation of
#' the statistic charted. This is based on results from Burr (1969).
#'
#' The function `limits.xbar` returns a matrix with lower and upper
#' control limits.
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @references Burr, I.W. (1969) Control charts for measurements with varying
#' sample sizes. *Journal of Quality Technology*, 1(3), 163-167.
#'
#' Montgomery, D.C. (2013) *Introduction to Statistical Quality Control*,
#' 7th ed. New York: John Wiley & Sons.
#'
#' Wetherill, G.B. and Brown, D.W. (1991) *Statistical Process Control*.
#' New York: Chapman & Hall.
#' @keywords htest hplot
stats.xbar <- function(data, sizes)
{
  data <- as.matrix(data)
  if(missing(sizes))
    sizes <- as.integer(rowSums(!is.na(data)))
  statistics <- rowMeans(data, na.rm = TRUE)
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.xbar <- function(data, sizes, std.dev = c("UWAVE-R", "UWAVE-SD", "MVLUE-R", "MVLUE-SD", "RMSDF"), ...)
{
  data <- as.matrix(data)
  if(missing(sizes))
    sizes <- as.integer(rowSums(!is.na(data)))
  if(any(sizes == 1))
    stop("group sizes must be larger than one")
  if(!is.numeric(std.dev))
    std.dev <- match.arg(std.dev, choices = eval(formals(sd.xbar)$std.dev))
  if(is.numeric(std.dev))
    { sd <- std.dev }
  else
    { switch(std.dev, 
             "UWAVE-R" = {  R <- apply(data, 1, function(x) 
                                       diff(range(x, na.rm = TRUE)))
                            d2 <- qcc.options("exp.R.unscaled")[sizes]
                            sd <- sum(R/d2)/length(sizes) 
                         }, 
             "UWAVE-SD" = { S <- apply(data, 1, sd, na.rm = TRUE)
                            sd <- sum(S/qcc.c4(sizes))/length(sizes) 
                          },
             "MVLUE-R"  = { R <- apply(data, 1, function(x) 
                            diff(range(x, na.rm = TRUE)))
                            d2 <- qcc.options("exp.R.unscaled")[sizes]
                            d3 <- qcc.options("se.R.unscaled")[sizes]
                            w  <- (d2/d3)^2
                            sd <- sum(R/d2*w)/sum(w) 
                          }, 
             "MVLUE-SD" = { S <- apply(data, 1, sd, na.rm = TRUE)
                            w  <- qcc.c4(sizes)^2/(1-qcc.c4(sizes)^2)
                            sd <- sum(S/qcc.c4(sizes)*w)/sum(w) 
                          },
             "RMSDF" =    { S <- apply(data, 1, sd, na.rm = TRUE)
                            w  <- sizes-1
                            sd <- sqrt(sum(S^2*w)/sum(w))/qcc.c4(sum(w)+1) 
                          }
      )
    }
  return(sd)
}

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
  .construct_limits(lcl,ucl)
}
