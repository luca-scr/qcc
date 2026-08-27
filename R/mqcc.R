#' Multivariate Quality Control Charts
#'
#' Create an object of class `'mqcc'` to perform multivariate statistical
#' quality control.
#'
#'
#' @param data For subgrouped data, a list with a data frame or a matrix for
#'   each variable to monitor. Each row of the data frame or matrix refers to a
#'   sample or ''rationale'' group.  For individual observations, where each
#'   sample has a single observation, users can provide a list with a data frame
#'   or a matrix having a single column, or a data frame or a matrix where each
#'   rows refer to samples and columns to variables. See examples.
#' @param type a character string specifying the type of chart:
#'   - `"T2"`: Hotelling \eqn{T^2} chart for subgrouped data.
#'   - `"T2.single"`: Hotelling \eqn{T^2} chart for individual observations.
#' @param center a vector of values to use for center of input variables.
#' @param cov a matrix of values to use for the covariance matrix of input
#'   variables.
#' @param limits a logical indicating if control limits (Phase I) must be
#'   computed (by default using [limits.T2()] or
#'   [limits.T2.single()]) and plotted, or a two-values vector
#'   specifying control limits.
#' @param pred.limits a logical indicating if prediction limits (Phase II) must
#'   be computed (by default using [limits.T2()] or
#'   [limits.T2.single()]) and plotted, or a two-values vector
#'   specifying prediction limits.
#' @param data.name a string specifying the name of the variable which appears
#'   on the plots. If not provided is taken from the object given as data.
#' @param labels a character vector of labels for each group.
#' @param newdata a data frame, matrix or vector, as for the `data`
#'   argument, providing further data to plot but not included in the
#'   computations.
#' @param newlabels a character vector of labels for each new group defined in
#'   the argument `newdata`.
#' @param confidence.level a numeric value between 0 and 1 specifying the
#'   confidence level of the computed probability limits.  By default is set at
#'   \eqn{(1 - 0.0027)^p} where \eqn{p} is the number of variables, and
#'   \eqn{0.0027} is the probability of Type I error for a single Shewhart chart
#'   at the usual 3-sigma control level.
#' @param plot logical. If `TRUE` a quality chart is plotted.
#' See `help(par)`.
#' @param digits the number of significant digits to use when `add.stats =
#' TRUE`.
#' @param restore.par a logical value indicating whether the previous
#' `par` settings must be restored. If you need to add points, lines, etc.
#' to a control chart set this to `FALSE`.
#' @param x an object of class `'mqcc'`.
#' @param object an object of class `'mqcc'`.
#' @param ... additional arguments to be passed to the generic function.
#' @return Returns an object of class `'mqcc'`.
#' @author Luca Scrucca
#' @seealso [stats.T2()], [stats.T2.single()], [limits.T2()], [limits.T2.single()], [ellipseChart()], [qcc()]
#' @references `r refs("mason_young_2002", "montgomery2013", "ryan_2011", "scrucca_2004", "wetherill_brown_1991")`
#' @export
#' @examples
#' ##  Subgrouped data
#' q  = mqcc(RyanMultivar, type = "T2")
#' summary(q)
#' ellipseChart(q)
#' ellipseChart(q, show.id = TRUE)
#' q  = mqcc(RyanMultivar, type = "T2", pred.limits = TRUE)
#'
#' # Xbar-charts for single variables computed adjusting the 
#' # confidence level of the T^2 chart:
#' q1  = with(RyanMultivar, 
#'            qcc(X1, type = "xbar", confidence.level = q$confidence.level^(1/2)))
#' summary(q1)
#' q2  = with(RyanMultivar,
#'            qcc(X2, type = "xbar", confidence.level = q$confidence.level^(1/2)))
#' summary(q2)
#'
#' require(MASS)
#' # generate new "in control" data
#' Xnew  = list(X1 = matrix(NA, 10, 4), X2 =  matrix(NA, 10, 4))
#' for(i in 1:4)
#'    { x  = mvrnorm(10, mu = q$center, Sigma = q$cov)
#'      Xnew$X1[,i]  = x[,1]
#'      Xnew$X2[,i]  = x[,2] 
#'    }
#' qq  = mqcc(RyanMultivar, type = "T2", newdata = Xnew, pred.limits = TRUE)
#' summary(qq)
#'
#' # generate new "out of control" data
#' Xnew  = list(X1 = matrix(NA, 10, 4), X2 =  matrix(NA, 10, 4))
#' for(i in 1:4)
#'    { x  = mvrnorm(10, mu = 1.2*q$center, Sigma = q$cov)
#'      Xnew$X1[,i]  = x[,1]
#'      Xnew$X2[,i]  = x[,2] 
#'    }
#' qq  = mqcc(RyanMultivar, type = "T2", newdata = Xnew, pred.limits = TRUE)
#' summary(qq)
#'
#' ## Individual observations data
#' q = mqcc(boiler, type = "T2.single", confidence.level = 0.999)
#' summary(q)
#'
#' # generate new "in control" data
#' boilerNew  = mvrnorm(10, mu = q$center, Sigma = q$cov)
#' qq  = mqcc(boiler, type = "T2.single", confidence.level = 0.999, 
#'            newdata = boilerNew, pred.limits = TRUE)
#' summary(qq)
#'
#' # generate new "out of control" data
#' boilerNew  = mvrnorm(10, mu = 1.01*q$center, Sigma = q$cov)
#' qq  = mqcc(boiler, type = "T2.single", confidence.level = 0.999, 
#'            newdata = boilerNew, pred.limits = TRUE)
#' summary(qq)
#'
#' # provides "robust" estimates of means and covariance matrix
#' rob  = cov.rob(boiler)
#' qrob  = mqcc(boiler, type = "T2.single", center = rob$center, cov = rob$cov)
#' summary(qrob)
mqcc <- function(data, type = c("T2", "T2.single"), center, cov,
                 limits = TRUE, pred.limits = FALSE,
                 data.name, labels, newdata, newlabels, 
                 confidence.level = (1-0.0027)^p, 
                 plot = TRUE, ...)
{
  call <- match.call()
  type <- match.arg(type)
  if(missing(data))
     stop("'data' argument is not specified")
  if(missing(data.name))
     data.name <- deparse(substitute(data))
  if(is.matrix(data))
     data <- as.data.frame(data)
  if(is.data.frame(data) | is.list(data))
       { data <- lapply(data, data.matrix) }
  else { stop("invalid data type!") }
  m <- unique(sapply(data, nrow))    # num. of samples
  if(length(m) > 1)
     stop("varying number of samples (rows)")
  n <- unique(sapply(data, ncol))   # samples sizes
  if(length(n) > 1)
     stop("varying sample size (columns)")
  if(n == 1) 
     type <- "T2.single"
  p <- length(data)                 # num. of variables
  var.names <- names(data)
  if(is.null(var.names))
    { var.names <- paste(data.name, "[", 1:p, "]", sep="") }
  #
  if(confidence.level <= 0 | confidence.level >= 1)
     stop("confidence.level must be a numeric value in the range (0,1)")
  if(missing(labels))
    { labels <- unique(unlist(sapply(data, rownames)))
      if(is.null(labels)) labels <- 1:m
      if(length(labels) != m)
         stop("labels must match the length of samples provided") }
  #
  if(missing(center)) center <- NULL
  if(missing(cov))    cov <- NULL
  #
  stats <- paste("stats.", type, sep = "")
  if(!exists(stats, mode="function"))
     stop(paste("function", stats, "is not defined"))
  stats <- do.call(stats, list(data, center = center, cov = cov))
  statistics <- stats$statistics
  stopifnot(length(labels) == length(statistics))
  names(statistics) <-  labels
  names(stats$center) <- var.names
  dimnames(stats$cov) <- list(var.names, var.names)

  # create object of class 'mqcc'
  object <- list(call = call, data.name = data.name,
                 var.names = var.names, data = data,
                 type = type,  confidence.level = confidence.level,
                 statistics = statistics,  means = stats$means,
                 center = stats$center, cov = stats$cov)
  class(object) <- "mqcc"

  # check for new data provided and update object
  if(!missing(newdata))
    { 
      newdata.name <- deparse(substitute(newdata))
      if(is.matrix(newdata))
         newdata <- as.data.frame(newdata)
      if(is.data.frame(newdata) | is.list(newdata))
           { newdata <- lapply(newdata, data.matrix) }
      else { stop("invalid data type!") }
      if(length(newdata) != p)
         stop("num. of variables for newdata not equal to num. of variables for data")
      stats <- paste("stats.", type, sep = "")
      if(!exists(stats, mode="function"))
         stop(paste("function", stats, "is not defined"))
      newstats <- do.call(stats, list(newdata, 
                                      center = object$center, 
                                      cov = object$cov))
      if(missing(newlabels))
        { start <- length(statistics)
          newlabels <- seq(start+1, start+length(newstats$statistics))
        }
      if(length(newlabels) != length(newstats$statistics))
        stop("labels must match the length of samples provided") 
      object$newdata  <- newdata
      object$newdata.name <- newdata.name
      names(newstats$statistics) <- newlabels
      object$newstats <- newstats$statistics
      object$newmeans <- newstats$means
      statistics <- c(statistics, newstats)
    }
  # compute control limits
  if(is.logical(limits))
    { if(limits)
        { limits <- paste("limits.", type, sep = "")
          if(!exists(limits, mode="function"))
            stop(paste("function", limits, "is not defined"))
          limits <- do.call(limits, list(ngroups = m, size = n, nvars = p, 
                                         conf = confidence.level))$control }
      else limits <- NULL                                    
    }
  else 
    { if(!is.numeric(limits))
         stop("'limits' must be a vector of length 2 or a 2-columns matrix")
      limits <- new_limits(matrix(limits, ncol = 2))
    }
  object$limits <- limits
  # compute prediction limits
  if(is.logical(pred.limits))
    { if(pred.limits)
        { pred.limits <- paste("limits.", type, sep = "")
          if(!exists(pred.limits, mode="function"))
            stop(paste("function", pred.limits, "is not defined"))
          pred.limits <- do.call(pred.limits, list(ngroups = m, size = n, nvars = p, 
                                                   conf = confidence.level))$prediction }
      else pred.limits <- NULL                                    
    }
  else 
    { if(!is.numeric(pred.limits))
         stop("'pred.limits' must be a vector of length 2 or a 2-columns matrix")
      pred.limits <- matrix(pred.limits, ncol = 2) |> new_limits(names = c("LPL", "UPL"))
    }
  object$pred.limits  <- pred.limits

  # identify violating rules observations
  violations <- rep(NA, length(statistics))
  if(is.numeric(object$limits))
  {
    violations[qccRulesViolatingWER1(object, limits = object$limits)] <- 1
  }
  if(is.numeric(object$pred.limits))
  {
    violations[qccRulesViolatingWER1(object, limits = object$pred.limits)] <- 2
  }
  object$violations <- violations
    
  if(plot) plot(object, ...) 
  #
  return(object)
}

#' @rdname mqcc
#' @method print mqcc
#' @export
#' @export print.mqcc
print.mqcc <- function(x, digits = getOption("digits"), ...)
{
  object <- x  # Argh.  Really want to use 'object' anyway
  # cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  cat(cli::rule(left = crayon::bold("Multivariate Quality Control Chart"), 
                width = min(getOption("width"),50)), "\n\n")

  data <- object$data
  m <- unique(sapply(data, nrow))     # num. of samples
  n <- unique(sapply(data, ncol))     # samples sizes
  p <- length(data)                   # num. of variables
  data.name <- object$data.name
  type <- object$type
  
  cat("Chart type                 =", type, "\n")
  cat("Data (phase I)             =", data.name, "\n")
  cat("Number of groups           =", m, "\n")
  cat("Group sample size          =", n, "\n")
  cat("Center = \n")
  print(object$center, digits = digits)
  cat("Covariance matrix = \n")
  print(object$cov, digits = digits)
  cat("|S| =", format(det(object$cov), digits = digits), "\n")
  
  newdata.name <- object$newdata.name
  newstats <- object$newstats
  if(!is.null(newstats)) 
  { 
    newdata <- object$newdata
    m <- unique(sapply(newdata, nrow))     # num. of samples
    n <- unique(sapply(newdata, ncol))     # samples sizes
    # p <- length(newdata)                   # num. of variables
    cat("\n")
    cat("New data (phase II)        =", newdata.name, "\n")
    cat("Number of groups           =", m, "\n")
    cat("Group sample size          =", n, "\n")
  }
  
  ctrl.limits <- object$limits
  if(!is.null(ctrl.limits)) 
  { 
    cat("\nControl limits:\n")
    .printShortMatrix(ctrl.limits, digits = digits, ...) 
  }

  pred.limits <- object$pred.limits
  if(!is.null(pred.limits)) 
  { 
    cat("\nPrediction limits:\n")
    .printShortMatrix(pred.limits, digits = digits, ...) 
  }
  
  invisible()        
}

#' @rdname mqcc
#' @method summary mqcc
#' @export
#' @export summary.mqcc
summary.mqcc <- function(object, ...) print.mqcc(object, ...)

#' @rdname mqcc
#' @method plot mqcc
#' @export
#' @export plot.mqcc
#' @inheritParams plot_common add.stats chart.all fill xlab ylab ylim
#' @param title a character string specifying the main title. Set `title =
#'   FALSE` or `title = NA` to remove the title.
#' @param label.limits a character vector specifying the labels for control
#'   limits (Phase I).
#' @param label.pred.limits a character vector specifying the labels for
#'   prediction control limits (Phase II).
#' @param axes.las numeric in \{0,1,2,3\} specifying the style of axis labels.
plot.mqcc <- function(x, 
                      add.stats = getOption("qcc.add.stats"),
                      chart.all = getOption("qcc.chart.all"),
                      fill = getOption("qcc.fill"),
                      label.limits = c("LCL", "UCL"),
                      label.pred.limits = c("LPL", "UPL"),
                      title, xlab, ylab, ylim, axes.las = 0,
                      digits = getOption("digits"),
                      restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  rules <- getOption("qcc.rules")
  zones <- getOption("qcc.zones")
  bg.margin <- getOption("qcc.bg.margin")
  bg.figure <- getOption("qcc.bg.figure")
  cex <- getOption("qcc.cex")
  font.stats <- getOption("qcc.font.stats")
  cex.stats <- getOption("qcc.cex.stats")
  if((missing(object)) | (!inherits(object, "mqcc")))
     stop("an object of class `mqcc' is required")
  # collect info from object
  data <- object$data
  m <- unique(sapply(data, nrow))     # num. of samples
  n <- unique(sapply(data, ncol))     # samples sizes
  p <- length(data)                   # num. of variables
  type <- object$type
  data.name <- object$data.name
  var  <- det(object$cov)
  stats <- object$statistics
  limits <- object$limits 
  pred.limits <- object$pred.limits 
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
  #
  if(chart.all) 
    { statistics <- c(stats, newstats)
      indices <- 1:length(statistics) }
  else
     { if(is.null(newstats))
         { statistics <- stats
           indices <- 1:length(statistics) }
       else
         { statistics <- newstats 
           indices <- seq(length(stats)+1, length(stats)+length(newstats)) }
     }
  #
  if(missing(title))
  { 
    if(is.null(newstats))
         title <- paste(type, "chart for", data.name)
       else if(chart.all)
               title <- paste(type, "chart for", data.name, "and", newdata.name)
            else title <- paste(type, "chart for", newdata.name) 
  }
  if(isFALSE(title) | is.na(title)) title <- ""
  
  cex.labels <- par("cex")*cex
  cex.stats.lines <- par("cex")*cex.stats
  
  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  par(bg  = bg.margin,
      cex = oldpar$cex * cex,
      mar = pmax(par("mar"), c(4.1,4.1,1.1,2.1), na.rm=TRUE),
      oma = if(add.stats) c(3.5*cex.stats.lines, 0, 1.5*cex.labels, 0)
            else          c(0, 0, 1.5*cex.labels, 0))

  # plot Shewhart chart
  plot(indices, statistics, type = "n",
       ylim = if(!missing(ylim)) ylim 
              else range(statistics, limits, pred.limits),
       ylab = if(missing(ylab)) "Group summary statistics" else ylab,
       xlab = if(missing(xlab)) "Group" else xlab, 
       axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = bg.figure)
  axis(1, at = indices, las = axes.las,
       labels = names(statistics) %||% as.character(indices),
       cex.axis = par("cex.axis")*0.9)
  axis(2, las = axes.las, cex.axis = par("cex.axis")*0.9)
  box()
  mtext(title, side = 3, outer = TRUE, 
        line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = par("cex")*cex,
        col  = par("col.main"))

  # draw prediction or control limits
  if(!is.null(pred.limits))
  {
    if(fill)
    { # fill the in-control area
      polygon(par("usr")[c(1,2,2,1)], pred.limits[c(1,1,2,2)], border = FALSE, 
              col = adjustcolor(zones$fill, alpha.f = 0.1))
    } else
    { 
      abline(h = pred.limits,
             lty = zones$lty[1],
             col = zones$col[1])
    }
    mtext(label.pred.limits, side = 4, at = c(pred.limits[1], pred.limits[2]), 
          las = 1, line = 0.1, col = gray(0.3),
          cex = par("cex")*cex.stats)
  }
  if(!is.null(limits))
  {
    if(fill)
    { # fill the in-control area
      polygon(par("usr")[c(1,2,2,1)], limits[c(1,1,2,2)], border = FALSE, 
              col = adjustcolor(zones$fill, alpha.f = 0.1))
    } else
    { 
      abline(h = limits,
             lty = zones$lty[1],
             col = zones$col[1])
    }
    mtext(label.limits, side = 4, at = c(limits[1], limits[2]), 
          las = 1, line = 0.1, col = gray(0.3),
          cex = par("cex")*cex.stats)
  }
  
  # draw lines & points
  lines(indices, statistics, type = "b", pch = NA)
  col <- rep(palette()[1], length(indices))
  pch <- rep(20, length(indices))
  if(!is.null(violations))
  { 
    i <- indices %in% which(violations > 0)
    col[i] <- rules$col[1]
    pch[i] <- rules$pch[1]
  }
  points(indices, statistics, col = col, pch = pch)

  if(chart.all & (!is.null(newstats)))
  { 
    len.obj.stats <- length(object$statistics)
    len.new.stats <- length(statistics) - len.obj.stats
    abline(v = len.obj.stats + 0.5, lty = 3)
    mtext("Calibration data", cex = par("cex")*0.8,
          at = len.obj.stats/2, line = 0, adj = 0.5)
    mtext("New data", cex = par("cex")*0.8, 
          at = len.obj.stats + len.new.stats/2, line = 0, adj = 0.5)
  }

  if(add.stats) 
  { 
    at <- c(0.10,0.40,0.65) 
    mtext(paste("Number of groups = ", length(statistics), sep = ""), 
          side = 1, outer = TRUE, line = 0*cex.stats.lines, adj = 0, at = at[1],
          font = font.stats,
          cex = par("cex")*cex.stats)
    mtext(paste("Sample size = ", signif(n, digits), sep = ""),
          side = 1, outer = TRUE, line = 1*cex.stats.lines, adj = 0, at = at[1],
          font = font.stats,
          cex = par("cex")*cex.stats)
    mtext(paste("|S| = ", signif(var, digits), sep = ""),
          side = 1, outer = TRUE, line = 2*cex.stats.lines, adj = 0, at = at[1],
          font = font.stats,
          cex = par("cex")*cex.stats)
    #
    if(is.numeric(limits))
    { 
      mtext(paste(label.limits[1], " = ", signif(limits[1], digits), sep = ""), 
            side = 1, outer = TRUE, line = 0*cex.stats.lines, adj = 0, at = at[2],
            font = font.stats,
            cex = par("cex")*cex.stats)
      mtext(paste(label.limits[2], " = ", signif(limits[2], digits), sep = ""),
            side = 1, outer = TRUE, line = 1*cex.stats.lines, adj = 0, at = at[2],
            font = font.stats,
            cex = par("cex")*cex.stats)
      mtext(paste("Num. beyond limits =",
                  sum(violations==1, na.rm=TRUE)), 
            side = 1, outer = TRUE, line = 2*cex.stats.lines, adj = 0, at = at[2],
            font = font.stats,
            cex = par("cex")*cex.stats)
    }        
    #
    if(is.numeric(pred.limits))
    { 
      mtext(paste(label.pred.limits[1], " = ", signif(pred.limits[1], digits), sep = ""),
            side = 1, outer = TRUE, line = 0*cex.stats.lines, adj = 0, at = at[3],
            font = font.stats,
            cex = par("cex")*cex.stats)
      mtext(paste(label.pred.limits[2], " = ", signif(pred.limits[2], digits), sep = ""),
            side = 1, outer = TRUE, line = 1*cex.stats.lines, adj = 0, at = at[3],
            font = font.stats,
            cex = par("cex")*cex.stats)
      mtext(paste("Num. beyond limits =",
                  sum(violations==2, na.rm=TRUE)), 
            side = 1, outer = TRUE, line = 2*cex.stats.lines, adj = 0, at = at[3],
            font = font.stats,
            cex = par("cex")*cex.stats)
    }        
  }
  
  invisible() 
}



#' Multivariate Quality Control Charts
#'
#' Plot an ellipse chart for a bivariate quality control data.
#'
#'
#' @inheritParams plot_common chart.all xlab ylab xlim ylim
#' @param object an object of class `'mqcc'`.
#' @param show.id a logical value indicating whether to plot point labels
#' (`TRUE`) or symbols (`FALSE`) for group means.
#' @param ngrid a value for the size of the grid over which the ellipse is
#' evaluated.
#' @param confidence.level a numeric value between 0 and 1 specifying the
#' confidence level of the computed probability limits.
#' @param correct.multiple a logical value indicating whether to correct or not
#' for multiple comparisons.
#' @param title a character string specifying the main title. Set `title =
#' FALSE` or `title = NA` to remove the title.
#' @param restore.par a logical value indicating whether the previous
#' `par` settings must be restored. If you need to add points, lines, etc.
#' to a control chart set this to `FALSE`.
#' @param ... additional arguments to be passed to the generic
#' [points()] function.
#' @author Luca Scrucca
#' @seealso [mqcc()], [stats.T2()], [stats.T2.single()]
#' @references `r refs("mason_young_2002", "montgomery2013", "ryan_2011")`
#' @export
#' @examples
#' # See examples in help(mqcc)
ellipseChart <- function(object,
                         chart.all = getOption("qcc.chart.all"),
                         show.id = FALSE, ngrid = 50,
                         confidence.level, correct.multiple = TRUE,
                         title, xlim, ylim, xlab, ylab,
                         restore.par = TRUE, ...) 
{
  rules <- getOption("qcc.rules")
  bg.margin <- getOption("qcc.bg.margin")
  bg.figure <- getOption("qcc.bg.figure")
  cex <- getOption("qcc.cex")
  if((missing(object)) | (!inherits(object, "mqcc")))
     stop("an object of class `mqcc' is required")

  data <- object$data
  m <- unique(sapply(data, nrow))     # num. of samples
  n <- unique(sapply(data, ncol))     # samples sizes
  p <- length(data)                   # num. of variables
  if(p > 2)
     stop("ellipse chart only available for bivariate data")
  center <- object$center
  cov <- object$cov
  # stats to plot: within-sample means   
  if(chart.all) 
    { stats <- rbind(object$means,object$newmeans) }
  else
    { if(is.null(object$newdata))
         stats <- object$means
      else
         stats <- object$newmeans }
  #
  if(missing(confidence.level))
     confidence.level <- object$confidence.level
  #
  alpha <- 1 - confidence.level
  if(correct.multiple) alpha <- 1-sqrt(1-alpha)
  # 
  # compute control limits for univariate charts
  if(all(n == 1))
    { q1 <- qcc(object$data[[1]], type="xbar.one", plot=FALSE,
                confidence.level = 1-alpha)
      q2 <- qcc(object$data[[2]], type="xbar.one", plot=FALSE,
                confidence.level = 1-alpha) 
  }
  else 
    { q1 <- qcc(object$data[[1]], type="xbar", plot=FALSE,
                confidence.level = 1-alpha)
      q2 <- qcc(object$data[[2]], type="xbar", plot=FALSE,
                confidence.level = 1-alpha) 
  }
  #
  if(missing(xlim))
     xlim <- range(pretty(stats[,1],1), q1$limits)
  if(missing(ylim))
     ylim <- range(pretty(stats[,2],1), q2$limits)
  if(missing(title))
  { 
    if(is.null(object$newstats))
      title <- paste("Ellipse chart for", object$data.name)
    else if(chart.all)
      title <- paste("Ellipse chart for", object$data.name, "and", object$newdata.name)
    else title <- paste("Ellipse chart for", object$newdata.name) 
  }
  if(isFALSE(title) | is.na(title)) title <- ""
  #
  grid <- cbind(seq(xlim[1], xlim[2], length = ngrid),
                seq(ylim[1], ylim[2], length = ngrid))
  grid.points <- expand.grid(grid[,1], grid[,2])
  T2 <- n * stats::mahalanobis(grid.points, center, cov)
  T2 <- matrix(T2, ngrid, ngrid)
  q <- object$limits[2]
  
  cex.labels <- par("cex")*cex
  
  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  par(bg  = bg.margin,
      cex = oldpar$cex * cex,
      mar = pmax(par("mar"), c(4.1,4.1,1.1,2.1), na.rm=TRUE),
      oma = c(0, 0, 1.5*cex.labels, 0))

  # plot ellipse chart
  plot(stats, type = "n", xlim = xlim, ylim = ylim, 
       ylab = if(missing(ylab)) object$var.names[2] else ylab,
       xlab = if(missing(xlab)) object$var.names[1] else xlab)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = bg.figure)
  box()
  mtext(title, side = 3, outer = TRUE, 
        line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = par("cex")*cex,
        col  = par("col.main"))
  
  contour(grid[,1], grid[,2], T2, levels = q, drawlabels = FALSE, add=TRUE)
  points(center[1], center[2], pch=3, cex=2)

  v <- which(object$violations > 0)
  col    <- rep(palette()[1], nrow(stats))
  col[v] <- rules$col[1]
  pch    <- rep(1, nrow(stats))
  pch[v] <- rules$pch[1]
  if(show.id) 
  { 
    text(stats, labels = names(object$statistics), 
         cex = 0.8*cex, col = col)
  } else        
  { 
    points(stats, col = col, pch = pch) 
  }

  if(!is.null(q1) & !is.null(q2)) 
  { 
    abline(v = q1$limits, lty = 2)
    abline(h = q2$limits, lty = 2) 
  }
  #
  invisible() 
}
