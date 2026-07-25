#' Statistics used in computing and drawing the Hotelling T^2 chart for individual observations data
#'
#' These functions are used to compute statistics required by the \eqn{T^2}
#' chart for individual observations.
#'
#' @export
#' @param data the observed data values
#' @param center a vector of values to use for center of input variables.
#' @param cov a matrix of values to use for the covariance matrix of input
#' variables.
#' @return The function `stats.T2.single` returns a list with components:
#' - `statistics`: a vector of values for the \eqn{T^2} statistic.
#' - `means`: a matrix of within group means for each variable, equal to
#'   `data` since samples are of size one.
#' - `center`: sample/group center statistic.
#' - `S`: covariance matrix.
#'
#' The function `limits.T2.single` returns a list with components:
#' - `control`: control limits.
#' - `prediction`: prediction limits.
#' @author Luca Scrucca
#' @seealso [mqcc()], [stats.T2()]
#' @references `r refs("mason_young_2002", "montgomery2013", "ryan_2011")`
stats.T2.single <- function(data, center = NULL, cov = NULL)
{ 
  data <- as.matrix(as.data.frame(data))
  m <- nrow(data)                       # num. of samples
  p <- ncol(data)                       # num. of variables
  n <- 1                                # samples sizes
  if(is.null(center))
    { center <- colMeans(data) }  # overall mean
  x <- scale(data, center = center, scale = FALSE)
  if(is.null(cov))
    { cov <- crossprod(x)/(m-1) }       # sample covar matrix
  cov.inv <- solve(cov)
  # Hotelling's T^2 statistic
  T2 <- apply(x, 1, function(x) x %*% cov.inv %*% x)
  list(statistics = T2, means = data, center = center, cov = cov)
}


#' @rdname stats.T2.single
#' @export
#' @param ngroups number of groups
#' @param size sample size
#' @param nvars number of variables
#' @param conf confidence level (0 < `conf` < 1)
limits.T2.single <- function(ngroups, size = 1, nvars, conf)
{ 
  m <- ngroups
  n <- size
  p <- nvars
  # Phase 1 control limits
  # Tracy Mason Young (1992)
  ucl <- (m-1)^2/m*qbeta(conf, p/2, (m-p-1)/2)
  # Phase 2 prediction limits
  upl <- p*(m+1)*(m-1)/(m*(m-p))*qf(conf, p, m-p)

  list(
    control = new_limits(0,ucl),
    prediction = new_limits(0,upl, names = c("LPL", "UPL"))
  )
}
