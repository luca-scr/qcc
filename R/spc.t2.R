#' Statistics used in computing and drawing the Hotelling T^2 chart for subgrouped data
#'
#' These functions are used to compute statistics required by the \eqn{T^2}
#' chart.
#'
#' @export
#' @param data the observed data values
#' @param center a vector of values to use for center of input variables.
#' @param cov a matrix of values to use for the covariance matrix of input
#' variables.
#' @return The function `stats.T2` returns a list with components:
#' - `statistics`: a vector of values for the \eqn{T^2} statistic.
#' - `means`: a matrix of within group means for each variable.
#' - `center`: sample/group center statistic.
#' - `S`: covariance matrix.
#'
#' The function `limits.T2` returns a list with components:
#' - `control`: control limits.
#' - `prediction`: prediction limits.
#' @author Luca Scrucca
#' @seealso [mqcc()], [stats.T2.single()]
#' @references `r refs("mason_young_2002", "montgomery2013", "ryan_2011")`
stats.T2 <- function(data, center = NULL, cov = NULL)
{ 
  data <- lapply(data, data.matrix) 
  m <- unique(sapply(data, nrow))[1]    # num. of samples
  n <- unique(sapply(data, ncol))       # samples sizes
  p <- length(data)                     # num. of variables
  #
  means <- lapply(data, function(x) 
                        rowMeans(x, na.rm = TRUE))  # within-sample means
  means <- as.matrix(as.data.frame(means))
  if(is.null(center))
     center <- sapply(data, mean, na.rm = TRUE)     # overall mean
  x <- scale(means, center = center, scale = FALSE)
  if(is.null(cov))
    { cov <- matrix(0, p, p)            # pooled within-sample covar matrix
      for(k in 1:m)
          cov <- cov + crossprod(scale(sapply(data, function(x) x[k,]),
                                       center = means[k,], scale = FALSE))/(n-1)
          # cov <- cov + var(sapply(data, function(x) x[k,]))
       cov <- cov/m 
    }
  cov.inv <- solve(cov)
  # Hotelling's T^2 statistic
  T2 <- n*apply(x, 1, function(x) x %*% cov.inv %*% x)
  list(statistics = T2, means = means, center = center, cov = cov)
}

#' @rdname stats.T2
#' @export
#' @param ngroups number of groups
#' @param size sample size
#' @param nvars number of variables
#' @param conf confidence level (0 < `conf` < 1)
limits.T2 <- function(ngroups, size, nvars,  conf)
{ 
  m   <- ngroups
  n   <- size
  p   <- nvars
  # Phase 1 control limits
  ucl <- p*(m-1)*(n-1)/(m*n-m-p+1)*qf(conf, p, m*n-m-p+1)
  lcl <- 0
  ctrl.limits <- matrix(c(lcl, ucl), ncol = 2)
  # Phase 2 prediction limits
  ucl <- p*(m+1)*(n-1)/(m*n-m-p+1)*qf(conf, p, m*n-m-p+1)
  lcl <- 0
  pred.limits <- matrix(c(lcl, ucl), ncol = 2)

  rownames(ctrl.limits) <- rownames(pred.limits) <- rep("", nrow(pred.limits))
  colnames(ctrl.limits) <- c("LCL", "UCL")
  colnames(pred.limits) <- c("LPL", "UPL")

  return(list(control = ctrl.limits, prediction = pred.limits))
}
