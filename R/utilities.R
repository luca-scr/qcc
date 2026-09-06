#' Grouping data based on a sample indicator
#'
#' This function allows to easily group data to use as input to the
#' `qcc()` function.
#'
#'
#' @param data a data frame (or a similar structure which can be coerced to be
#' a `data.frame`) providing the observed data. If not provided, the input
#' to the following arguments must be objects defined in the calling
#' environment.
#' @param x a name from `data` or a vector of observed data values.
#' @param sample a name from `data` or a vector of sample indicators
#' defining the rationale subgroups of data values.
#' @return The function returns a matrix of suitable dimensions. If one or more
#' group have fewer observations than others, `NA` values are used to fill
#' empty values.
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @export
#' @examples
#'
#' data(pistonrings)
#' # create a matrix of 40 samples made of 5 observations each
#' qccGroups(data = pistonrings, diameter, sample)
#' # remove some observations to get still a 40x5 matrix but filled with NAs 
#' qccGroups(data = pistonrings[-c(1,2,50,52,199),], diameter, sample)
#'
qccGroups <- function(data, x, sample)
{
  # collect x and sample from data if provided
  if(!missing(data))
  {
    x      <- eval(substitute(x), data, parent.frame())
    sample <- eval(substitute(sample), data, parent.frame())
  }
  stopifnot(length(x) == length(sample))
  #
  x <- lapply(split(x, sample), as.vector)
  lx <- sapply(x, length)
  for(i in which(lx != max(lx)))
      x[[i]] <- c(x[[i]], rep(NA, max(lx)-lx[i]))
  x <- t(sapply(x, as.vector))
  return(x)
}

#' Histogram number of classes/bins
#'
#' Computes the optimal number of classes/bins for an histogram as the maximum
#' between the Sturges and Freedman-Diaconis (FD) estimators. For small
#' datasets the Sturges value will usually be chosen, while larger datasets
#' will usually default to FD. Avoids the overly conservative behaviour of FD
#' and Sturges for small and large datasets, respectively. This is the default
#' option in `numpy.histogram_bin_edges` available in Python.
#'
#'
#' @param x a vector of data values.
#' @param ... additional arguments to be passed to low level functions.
#' @return The value of suggested number of classes/bins.
#' @author Luca Scrucca
#' @seealso [grDevices::nclass.FD()], [grDevices::nclass.Sturges()]
#' @export
#' @examples
#'
#' set.seed(1)
#' x <- stats::rnorm(111)
#' nclass.hist(x)
#' x <- stats::rnorm(1111)
#' nclass.hist(x)
#'
nclass.hist <- function(x, ...)
{
  max(nclass.FD(x), nclass.Sturges(x))
}



#' Overdispersion test for binomial and poisson data
#'
#' This function allows to test for overdispersed data in the binomial and
#' poisson case.
#'
#' This very simple test amounts to compute the test statistic
#' \deqn{D = s^2 / \sigma^2 \times (n - 1)}
#' where \eqn{s^2} is the observed variance, \eqn{\sigma^2} is the theoretical
#' variance, and \eqn{n} is the number of observations. The test statistic is
#' then compared to the critical value of a Chi-square distribution with
#' \eqn{n-1} degrees of freedom.
#'
#' @param x a vector of observed data values
#' @param size for binomial data, a vector of sample sizes
#' @param type a character string specifying the distribution for testing,
#' either `"poisson"` or `"binomial"`. By default, if `size` is
#' provided a binomial distribution is assumed, otherwise a poisson
#' distribution.
#' @return The function returns a matrix of results.
#' @author Luca Scrucca
#' @references `r refs("wetherill_brown_1991")`
#' @keywords htest
#' @export
#' @examples
#'
#' # data from Wetherill and Brown (1991) pp. 212--213, 216--218:
#' x  = c(12,11,18,11,10,16,9,11,14,15,11,9,10,13,12,
#'        8,12,13,10,12,13,16,12,18,16,10,16,10,12,14)
#' size  = rep(50, length(x))
#' qccOverdispersionTest(x, size)
#'
#' x  = c(11,8,13,11,13,17,25,23,11,16,9,15,10,16,12,
#'        8,9,15,4,12,12,12,15,17,14,17,12,12,7,16)
#' qccOverdispersionTest(x)
#'
qccOverdispersionTest <- function(x, size, 
                                  type = ifelse(missing(size), "poisson", "binomial"))
{
  type <- match.arg(type, c("poisson", "binomial"))
  if (type=="binomial" & missing(size))
     stop("binomial data require argument \"size\"")
  if (!missing(size))
     if (length(x) != length(size))   
        stop("arguments \"x\" and \"size\" must be vector of same length")

  n <- length(x)
  obs.var <- var(x)
  if (type=="binomial")
     { p <- sum(x)/sum(size)
       theor.var <- mean(size)*p*(1-p) }
  else if (type=="poisson")
          { theor.var <- mean(x) }
       else
          stop("invalid \"type\" argument. See help.")

  D <- (obs.var * (n-1)) / theor.var
  p.value <- 1-pchisq(D, n-1)

  out <- matrix(c(obs.var/theor.var, D, signif(p.value,5)), 1, 3)
  rownames(out) <- paste(type, "data")
  colnames(out) <- c("Obs.Var/Theor.Var", "Statistic", "p-value") 
  names(dimnames(out)) <- c(paste("Overdispersion test"), "")
  return(out)
}

#' Compute c4 constant
#'
#' Computes the c4 bias-correction factor for the sample standard deviation at
#' a given sample size.
#' @param n sample size(s)
#' @keywords internal
#' @noRd
# DEPRECATE qcc.c4 -> .c4
qcc.c4 <- function(n) sqrt(2/(n - 1)) * exp(lgamma(n/2) - lgamma((n - 1)/2))

#' Construct Control Limits
#'
#' Creates a two-column matrix of lower and upper control limits.
#'
#' @param ... A two-column matrix or separate lower and upper limit vectors.
#' @param names Column names for the lower and upper limits.
#'
#' @keywords internal
#' @noRd
new_limits <- function(..., names = c("LCL", "UCL")) {
  limits <- cbind(...)
  dimnames(limits) <- list(rep("", nrow(limits)), names)
  limits
}

#' Blue color palette
#'
#' Creates a palette of blue colors.
#'
#' @param n the number of colors.
#' @return A character vector containing the requested colors.
#' @export
blues.colors <- function (n) 
{
  palette <- grDevices::colorRampPalette(c("#03396c", "#005b96", "#6497b1", "#b3cde0"), 
                                         space = "Lab")
  palette(n)
}

#' Print Representative Columns and Rows for a Matrix
#'
#' Prints a shortened matrix containing configurable numbers of leading and
#' trailing rows and columns.
#'
#' @param x A matrix to summarize.
#' @param head number of initial rows to print
#' @param tail number of last rows to print
#' @param chead number of initial columns to print
#' @param ctail number of last columns to print
#' @param ... extra arguments passed to [print()]
#' @keywords internal
#' @noRd
.printShortMatrix <- function(x, head = 2, tail = 1, chead = 5, ctail = 1, ...)
{ 
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  rnames <- rownames(x)
  cnames <- colnames(x)
  dnames <- names(dimnames(x))
  
  if(is.na(head <- as.numeric(head))) head <- 2
  if(is.na(tail <- as.numeric(tail))) tail <- 1
  if(is.na(chead <- as.numeric(chead))) chead <- 5
  if(is.na(ctail <- as.numeric(ctail))) ctail <- 1
  
  if(nr > (head + tail))
  { 
    if(is.null(rnames)) 
      rnames <- paste("[", 1:nr, ",]", sep ="")
    x <- rbind(x[1:head,,drop=FALSE], 
               rep(NA, nc), 
               x[(nr-tail+1):nr,,drop=FALSE])
    rownames(x) <- c(rnames[1:head], ":", rnames[(nr-tail+1):nr])
  }
  if(nc > (chead + ctail))
  { 
    if(is.null(cnames)) 
      cnames <- paste("[,", 1:nc, "]", sep ="")
    x <- cbind(x[,1:chead,drop=FALSE], 
               rep(NA, nrow(x)), 
               x[,(nc-ctail+1):nc,drop=FALSE])
    colnames(x) <- c(cnames[1:chead], "...", cnames[(nc-ctail+1):nc])
  }
  names(dimnames(x)) <- dnames
  print(x, na.print = "", ...)
  invisible(x)
}

#' qcc Options Interface
#'
#' `qcc.options()` is retained as a deprecated compatibility interface for the
#' standard `qcc.*` options.
#'
#' Set plotting defaults with [options()] using these names:
#' `qcc.add.stats`, `qcc.chart.all`, `qcc.fill`, `qcc.rules`, `qcc.zones`,
#' `qcc.bg.margin`, `qcc.bg.figure`, `qcc.cex`, `qcc.font.stats`, and
#' `qcc.cex.stats`.
#'
#' @param ... No arguments to view all `qcc.*` options; a single legacy option
#'   name to retrieve its value; or named values or a named list to update
#'   standard options.
#' @return A named list for a call without arguments, a single option value for
#'   a character-name lookup, or the updated named list invisibly after a
#'   setter call.
#' @export
qcc.options <- function(...) {
  warning(
    paste0(
      "`qcc.options()` is deprecated; use `options(qcc.<name> = value)`."
    ),
    call. = FALSE
  )

  current_options <- function() {
    current <- options()
    current <- current[startsWith(names(current), "qcc.")]
    names(current) <- substring(names(current), 5L)
    current
  }

  current <- current_options()
  if(nargs() == 0L)
    return(current)

  values <- list(...)
  if(length(values) == 1L && is.null(names(values))) {
    value <- values[[1L]]
    switch(
      mode(value),
      list = values <- value,
      character = return(getOption(paste0("qcc.", value))),
      stop(paste("invalid argument:", sQuote(value)))
    )
  }
  if(length(values) == 0L)
    return(current)

  option_names <- names(values)
  if(is.null(option_names) || any(!nzchar(option_names)))
    stop("options must be given by name")

  names(values) <- paste0("qcc.", option_names)
  options(values)
  invisible(current_options())
}


#' Validate Sample Sizes
#'
#' Replace invalid sample sizes with NA, or error in strict mode.
#'
#' @keywords internal
#' @noRd
assert_n <- function(n, strict = FALSE) {
  invalid <- !is.finite(n) | n < 2 | n != floor(n)

  if (!any(invalid)) return(invisible(n))

  if (strict) cli::cli_abort("Invalid sample sizes: {n[invalid]}.")

  cli_warn("Replacing invalid sample sizes with `NA`: {n[invalid]}.")
  invisible(replace(n, invalid, NA_real_))
}
