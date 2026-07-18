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
#' @references Wetherill, G.B. and Brown, D.W. (1991) *Statistical Process
#' Control*, New York, Chapman and Hall, pp. 216--218
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
qcc.c4 <- \(n) sqrt(2/(n - 1)) * exp(lgamma(n/2) - lgamma((n - 1)/2))


#' Control Limits Constructor
#'
#' Returns lower and upper control limit vectors in a consistent structure.
#' Used by limits.*() functions.
#'
#' @keywords internal
.construct_limits <- function(lcl,ucl) {
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# FIX: blues.colors is used by ocCurves and should not be documented with paretoChart

#' @rdname paretoChart
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

#' Set or return options for the qcc package.
#'
#' This function can be used to control the behavior of the 'qcc' library such
#' as the background color, out-of-control points appearance, and many others.
#'
#' The available options are:
#'
#' - `exp.R.unscaled`: a vector specifying, for each sample size, the expected
#'   value of the relative range (i.e. \eqn{R/\sigma}) for a normal
#'   distribution. This appears as \eqn{d_2} on most tables containing factors
#'   for the construction of control charts.
#' - `se.R.unscaled`: a vector specifying, for each sample size, the standard
#'   error of the relative range (i.e. \eqn{R/\sigma}) for a normal
#'   distribution. This appears as \eqn{d_3} on most tables containing factors
#'   for the construction of control charts.
#' - `beyond.limits$pch`: plotting character used to highlight points beyond
#'   control limits.
#' - `beyond.limits$col`: color used to highlight points beyond control
#'   limits.
#' - `violating.runs$pch`: plotting character used to highlight points
#'   violating runs.
#' - `violating.runs$col`: color used to highlight points violating runs.
#' - `run.length`: the maximum value of a run before to signal a point as out
#'   of control.
#' - `bg.margin`: background color used to draw the margin of the charts.
#' - `bg.figure`: background color used to draw the figure of the charts.
#' - `cex`: character expansion used to draw plot annotations (labels, title,
#'   tickmarks, etc.).
#' - `font.stats`: font used to draw text at the bottom of control charts.
#' - `cex.stats`: character expansion used to draw text at the bottom of
#'   control charts.
#'
#' @param ... the option to be set or retrieved. See details.
#' @return If the functions is called with no argument return a list of
#' available options.
#'
#' If an option argument is provided the corresponding value is returned.
#'
#' If a value is associated with an option argument, such option is set and the
#' list of updated option values is invisibly returned. In this case the list
#' `.qcc.options` is modified and any modification will remain in effect
#' for the rest of the session.
#' @author Luca Scrucca
#' @seealso [qcc()]
#' @keywords htest hplot
#' @export
#' @examples
#'
#' old  = qcc.options()			# save defaults
#' qcc.options("cex.stats")		# get a single parameter
#' qcc.options("cex.stats"=1.2)	# change parameters
#' qcc.options(bg.margin="azure2")
#' qcc.options("violating.runs" = list(pch = 15, col = "purple"))
#' qcc.options("beyond.limits" = list(pch = 15, col = "orangered"))
#' qcc(rnorm(100), type = "xbar.one", std.dev = 0.7)	# see the results
#' qcc.options(old)				# restore old defaults 
#'
qcc.options <- function(...)
{
  current <- .qcc.options
  if(nargs() == 0) return(current)
#  if(is.character(...))
#       temp <- eval(parse(text = paste(c("list(", ..., ")"))))
#  else temp <- list(...)
  temp <- list(...)
  if(length(temp) == 1 && is.null(names(temp))) 
    { arg <- temp[[1]]
      switch(mode(arg),
             list = temp <- arg,
             character = return(.qcc.options[[arg]]),
             stop(paste("invalid argument:", sQuote(arg)))) }
  if(length(temp) == 0) return(current)
  name <- names(temp)
  if(is.null(name)) stop("options must be given by name")
  changed <- current[name]
  current[name] <- temp
  env <- if(sys.parent() == 0) asNamespace("qcc") 
         else                  parent.frame()
  assign(".qcc.options", current, envir = env)
  invisible(current)
}


#' Default `qcc` Settings
#'
#' Stores the statistical constants and graphical settings used as the qcc
#' package defaults. See [qcc.options()].
#'
#' @keywords internal
".qcc.options" <- list(
  exp.R.unscaled = c(NA, 1.128, 1.693, 2.059, 2.326, 2.534, 2.704, 2.847, 2.970, 3.078, 3.173, 3.258, 3.336, 3.407, 3.472, 3.532, 3.588, 3.640, 3.689, 3.735, 3.778, 3.819, 3.858, 3.895, 3.931),
  se.R.unscaled = c(NA, 0.8525033, 0.8883697, 0.8798108, 0.8640855, 0.8480442, 0.8332108, 0.8198378, 0.8078413, 0.7970584, 0.7873230, 0.7784873, 0.7704257, 0.7630330, 0.7562217, 0.7499188, 0.7440627, 0.7386021, 0.7334929, 0.7286980, 0.7241851, 0.7199267, 0.7158987, 0.7120802, 0.7084528, 0.7050004, 0.7017086, 0.6985648, 0.6955576, 0.6926770, 0.6899137, 0.6872596, 0.6847074, 0.6822502, 0.6798821, 0.6775973, 0.6753910, 0.6732584, 0.6711952, 0.6691976, 0.6672619, 0.6653848, 0.6635632, 0.6617943, 0.6600754, 0.6584041, 0.6567780, 0.6551950, 0.6536532, 0.6521506),
  rules = list(col = c("#F03B20", "#EE7600", "#FD8D3C", "#CD5555", "#7B3294", "#008837", "#B0BC17", "#C51B7D"),
                   # c("#F03B20", "#FD8D3C", "#FEB24C", "#FED976"), 
               pch = c(19, 15, 17, 8, 19, 15, 17, 18)),
  zones = list(fill = "#5E81AC", # "#81a1c1",
               lty = c(2,2,2), 
               col = grey(c(0.1, 0.4, 0.7))),
  bg.margin = "#EFF0F2", # grey(0.915),
  bg.figure = "white",
  cex = 1,
  font.stats = 1,
  cex.stats = 0.9,
  add.stats = TRUE,
  chart.all = TRUE, 
  fill = TRUE)
  
