#-----------------------------------------------------------------------------#
#                                                                             #
#  MULTIVARIATE CONTROL CHARTS                                                #
#                                                                             #
#  Written by: Luca Scrucca                                                   #
#              Department of Statistics                                       #
#              University of Perugia, ITALY                                   #
#              luca@stat.unipg.it                                             #
#                                                                             #
# Last modified: October 2009                                                 #
#-----------------------------------------------------------------------------#



#' Multivariate Quality Control Charts
#'
#' Create an object of class `'mqcc'` to perform multivariate statistical
#' quality control.
#'
#'
#' @aliases mqcc print.mqcc summary.mqcc plot.mqcc
#' @export mqcc
#' @param data For subgrouped data, a list with a data frame or a matrix for
#' each variable to monitor. Each row of the data frame or matrix refers to a
#' sample or ''rationale'' group.  For individual observations, where each
#' sample has a single observation, users can provide a list with a data frame
#' or a matrix having a single column, or a data frame or a matrix where each
#' rows refer to samples and columns to variables. See examples.
#' @param type a character string specifying the type of chart:
#'
#' | Type | Chart description |
#' | --- | --- |
#' | `"T2"` | Hotelling \eqn{T^2} chart for subgrouped data |
#' | `"T2.single"` | Hotelling \eqn{T^2} chart for individual observations |
#' @param center a vector of values to use for center of input variables.
#' @param cov a matrix of values to use for the covariance matrix of input
#' variables.
#' @param limits a logical indicating if control limits (Phase I) must be
#' computed (by default using [limits.T2()] or
#' [limits.T2.single()]) and plotted, or a two-values vector
#' specifying control limits.
#' @param pred.limits a logical indicating if prediction limits (Phase II) must
#' be computed (by default using [limits.T2()] or
#' [limits.T2.single()]) and plotted, or a two-values vector
#' specifying prediction limits.
#' @param data.name a string specifying the name of the variable which appears
#' on the plots. If not provided is taken from the object given as data.
#' @param labels a character vector of labels for each group.
#' @param newdata a data frame, matrix or vector, as for the `data`
#' argument, providing further data to plot but not included in the
#' computations.
#' @param newlabels a character vector of labels for each new group defined in
#' the argument `newdata`.
#' @param confidence.level a numeric value between 0 and 1 specifying the
#' confidence level of the computed probability limits.  By default is set at
#' \eqn{(1 - 0.0027)^p} where \eqn{p} is the number of variables, and
#' \eqn{0.0027} is the probability of Type I error for a single Shewhart chart
#' at the usual 3-sigma control level.
#' @param plot logical. If `TRUE` a quality chart is plotted.
#' @param add.stats a logical value indicating whether statistics and other
#' information should be printed at the bottom of the chart.
#' @param chart.all a logical value indicating whether both statistics for
#' `data` and for `newdata` (if given) should be plotted.
#' @param fill a logical value specifying if the in-control area should be
#' filled with the color specified in `qcc.options("zones")$fill`.
#' @param label.limits a character vector specifying the labels for control
#' limits (Phase I).
#' @param label.pred.limits a character vector specifying the labels for
#' prediction control limits (Phase II).
#' @param title a character string specifying the main title. Set `title =
#' FALSE` or `title = NA` to remove the title.
#' @param xlab a string giving the label for the x-axis.
#' @param ylab a string giving the label for the y-axis.
#' @param ylim a numeric vector specifying the limits for the y-axis.
#' @param axes.las numeric in {0,1,2,3} specifying the style of axis labels.
#' See `help(par)`.
#' @param digits the number of significant digits to use when `add.stats =
#' TRUE`.
#' @param restore.par a logical value indicating whether the previous
#' `par` settings must be restored. If you need to add points, lines, etc.
#' to a control chart set this to `FALSE`.
#' @param x an object of class `'mqcc'`.
#' @param ... additional arguments to be passed to the generic function.
#' @return Returns an object of class `'mqcc'`.
#' @author Luca Scrucca
#' @seealso [stats.T2()], [stats.T2.single()],
#' [limits.T2()], [limits.T2.single()],
#' [ellipseChart()], [qcc()]
#' @references Mason, R.L. and Young, J.C. (2002) *Multivariate
#' Statistical Process Control with Industrial Applications*, SIAM.
#'
#' Montgomery, D.C. (2013) *Introduction to Statistical Quality Control*,
#' 7th ed. New York: John Wiley & Sons.
#'
#' Ryan, T. P. (2011), *Statistical Methods for Quality Improvement*, 3rd
#' ed. New York: John Wiley & Sons, Inc.
#'
#' Scrucca, L. (2004). qcc: an R package for quality control charting and
#' statistical process control. *R News* 4/1, 11-17.
#'
#' Wetherill, G.B. and Brown, D.W. (1991) *Statistical Process Control*.
#' New York: Chapman & Hall.
#' @keywords htest hplot multivariate
#' @examples
#'
#' ##
#' ##  Subgrouped data
#' ##
#'
#' data(RyanMultivar)
#' str(RyanMultivar)
#'
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
#' ##
#' ## Individual observations data
#' ##
#'
#' data(boiler)
#' str(boiler)
#'
#' q  = mqcc(boiler, type = "T2.single", confidence.level = 0.999)
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
#'
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
      limits <- matrix(limits, ncol = 2)
      dimnames(limits) <- list(rep("",nrow(limits)), c("LCL ", "UCL"))
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
      pred.limits <- matrix(pred.limits, ncol = 2)
      dimnames(pred.limits) <- list(rep("",nrow(pred.limits)), c("LPL ", "UPL"))
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
plot.mqcc <- function(x, 
                      add.stats = qcc.options("add.stats"), 
                      chart.all = qcc.options("chart.all"), 
                      fill = qcc.options("fill"),
                      label.limits = c("LCL", "UCL"),
                      label.pred.limits = c("LPL", "UPL"),
                      title, xlab, ylab, ylim, axes.las = 0,
                      digits = getOption("digits"),
                      restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
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
  
  cex.labels <- par("cex")*qcc.options("cex")
  cex.stats <- par("cex")*qcc.options("cex.stats")
  
  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = pmax(par("mar"), c(4.1,4.1,1.1,2.1), na.rm=TRUE),
      oma = if(add.stats) c(3.5*cex.stats, 0, 1.5*cex.labels, 0) 
            else          c(0, 0, 1.5*cex.labels, 0))

  # plot Shewhart chart
  plot(indices, statistics, type = "n",
       ylim = if(!missing(ylim)) ylim 
              else range(statistics, limits, pred.limits),
       ylab = if(missing(ylab)) "Group summary statistics" else ylab,
       xlab = if(missing(xlab)) "Group" else xlab, 
       axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  axis(1, at = indices, las = axes.las,
       labels = names(statistics) %||% as.character(indices),
       cex.axis = par("cex.axis")*0.9)
  axis(2, las = axes.las, cex.axis = par("cex.axis")*0.9)
  box()
  mtext(title, side = 3, outer = TRUE, 
        line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = par("cex")*qcc.options("cex"), 
        col  = par("col.main"))

  # draw prediction or control limits
  if(!is.null(pred.limits))
  {
    if(fill)
    { # fill the in-control area
      polygon(par("usr")[c(1,2,2,1)], pred.limits[c(1,1,2,2)], border = FALSE, 
              col = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.1))
    } else
    { 
      abline(h = pred.limits,
             lty = qcc.options("zones")$lty[1], 
             col = qcc.options("zones")$col[1])
    }
    mtext(label.pred.limits, side = 4, at = c(pred.limits[1], pred.limits[2]), 
          las = 1, line = 0.1, col = gray(0.3),
          cex = par("cex")*qcc.options("cex.stats"))
  }
  if(!is.null(limits))
  {
    if(fill)
    { # fill the in-control area
      polygon(par("usr")[c(1,2,2,1)], limits[c(1,1,2,2)], border = FALSE, 
              col = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.1))
    } else
    { 
      abline(h = limits,
             lty = qcc.options("zones")$lty[1], 
             col = qcc.options("zones")$col[1])
    }
    mtext(label.limits, side = 4, at = c(limits[1], limits[2]), 
          las = 1, line = 0.1, col = gray(0.3),
          cex = par("cex")*qcc.options("cex.stats"))
  }
  
  # draw lines & points
  lines(indices, statistics, type = "b", pch = NA)
  col <- rep(palette()[1], length(indices))
  pch <- rep(20, length(indices))
  if(!is.null(violations))
  { 
    i <- indices %in% which(violations > 0)
    col[i] <- qcc.options("rules")$col[1]
    pch[i] <- qcc.options("rules")$pch[1]
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
    # TOREMOVE computes the x margins of the figure region
    # plt <- par()$plt; usr <- par()$usr
    # px <- diff(usr[1:2])/diff(plt[1:2])
    # xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
    # at.col <- xfig[1] + diff(xfig[1:2])*c(0.15, 0.45, 0.75)
    # top.line <- 4.5
    # write info at bottom
    mtext(paste("Number of groups = ", length(statistics), sep = ""), 
          side = 1, outer = TRUE, line = 0*cex.stats, adj = 0, at = at[1],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    mtext(paste("Sample size = ", signif(n, digits), sep = ""),
          side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[1],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    mtext(paste("|S| = ", signif(var, digits), sep = ""),
          side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[1],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    #
    if(is.numeric(limits))
    { 
      mtext(paste(label.limits[1], " = ", signif(limits[1], digits), sep = ""), 
            side = 1, outer = TRUE, line = 0*cex.stats, adj = 0, at = at[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste(label.limits[2], " = ", signif(limits[2], digits), sep = ""),
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Num. beyond limits =",
                  sum(violations==1, na.rm=TRUE)), 
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }        
    #
    if(is.numeric(pred.limits))
    { 
      mtext(paste(label.pred.limits[1], " = ", signif(pred.limits[1], digits), sep = ""),
            side = 1, outer = TRUE, line = 0*cex.stats, adj = 0, at = at[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste(label.pred.limits[2], " = ", signif(pred.limits[2], digits), sep = ""),
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Num. beyond limits =",
                  sum(violations==2, na.rm=TRUE)), 
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }        
  }
  
  invisible() 
}



#' Multivariate Quality Control Charts
#'
#' Plot an ellipse chart for a bivariate quality control data.
#'
#'
#' @param object an object of class `'mqcc'`.
#' @param chart.all a logical value indicating whether both statistics for
#' `data` and for `newdata` (if given) should be plotted.
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
#' @param xlim a numeric vector specifying the limits for the x-axis.
#' @param ylim a numeric vector specifying the limits for the y-axis.
#' @param xlab a string giving the label for the x-axis.
#' @param ylab a string giving the label for the y-axis.
#' @param restore.par a logical value indicating whether the previous
#' `par` settings must be restored. If you need to add points, lines, etc.
#' to a control chart set this to `FALSE`.
#' @param ... additional arguments to be passed to the generic
#' [points()] function.
#' @author Luca Scrucca
#' @seealso [mqcc()], [stats.T2()],
#' [stats.T2.single()]
#' @references Mason, R.L. and Young, J.C. (2002) *Multivariate
#' Statistical Process Control with Industrial Applications*, SIAM.
#'
#' Montgomery, D.C. (2013) *Introduction to Statistical Quality Control*,
#' 7th ed. New York: John Wiley & Sons.
#'
#' Ryan, T. P. (2011), *Statistical Methods for Quality Improvement*, 3rd
#' ed. New York: John Wiley & Sons, Inc.
#' @keywords htest hplot multivariate
#' @export
#' @examples
#'
#' # See examples in help(mqcc)
#'
ellipseChart <- function(object, chart.all = TRUE, show.id = FALSE, ngrid = 50,
                         confidence.level, correct.multiple = TRUE,
                         title, xlim, ylim, xlab, ylab,
                         restore.par = TRUE, ...) 
{
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
  x     <- grid - matrix(center, ngrid, 2, byrow=TRUE)
  cov.inv <- solve(object$cov)
  T2    <- n*apply(expand.grid(x[,1], x[,2]), 1,
                               function(x) x %*% cov.inv %*% x)
  T2    <- matrix(T2, ngrid, ngrid)
  q <- object$limits[2]
  
  cex.labels <- par("cex")*qcc.options("cex")
  cex.stats <- par("cex")*qcc.options("cex.stats")
  
  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = pmax(par("mar"), c(4.1,4.1,1.1,2.1), na.rm=TRUE),
      oma = c(0, 0, 1.5*cex.labels, 0))

  # plot ellipse chart
  plot(stats, type = "n", xlim = xlim, ylim = ylim, 
       ylab = if(missing(ylab)) object$var.names[2] else ylab,
       xlab = if(missing(xlab)) object$var.names[1] else xlab)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  box()
  mtext(title, side = 3, outer = TRUE, 
        line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = par("cex")*qcc.options("cex"), 
        col  = par("col.main"))
  
  contour(grid[,1], grid[,2], T2, levels = q, drawlabels = FALSE, add=TRUE)
  points(center[1], center[2], pch=3, cex=2)

  v <- which(object$violations > 0)
  col    <- rep(palette()[1], nrow(stats))
  col[v] <- qcc.options("rules")$col[1]
  pch    <- rep(1, nrow(stats))
  pch[v] <- qcc.options("rules")$pch[1]
  if(show.id) 
  { 
    text(stats, labels = names(object$statistics), 
         cex = 0.8*qcc.options("cex"), col = col)
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

# T2 chart



#' Statistics used in computing and drawing the Hotelling T^2 chart for
#' subgrouped data
#'
#' These functions are used to compute statistics required by the \eqn{T^2}
#' chart.
#'
#'
#' @aliases stats.T2 limits.T2
#' @export stats.T2
#' @export limits.T2
#' @param data the observed data values
#' @param center a vector of values to use for center of input variables.
#' @param cov a matrix of values to use for the covariance matrix of input
#' variables.
#' @param ngroups number of groups
#' @param size sample size
#' @param nvars number of variables
#' @param conf confidence level (0 < `conf` < 1)
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
#' @references Mason, R.L. and Young, J.C. (2002) *Multivariate
#' Statistical Process Control with Industrial Applications*, SIAM.
#'
#' Montgomery, D.C. (2013) *Introduction to Statistical Quality Control*,
#' 7th ed. New York: John Wiley & Sons.
#'
#' Ryan, T. P. (2011), *Statistical Methods for Quality Improvement*, 3rd
#' ed. New York: John Wiley & Sons, Inc.
#' @keywords htest hplot multivariate
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

limits.T2 <- function(ngroups, size, nvars,  conf)
{ 
  m   <- ngroups     # num. of samples
  n   <- size        # samples size
  p   <- nvars       # num. of variables
  # Phase 1 control limits   
  ucl <- p*(m-1)*(n-1)/(m*n-m-p+1)*qf(conf, p, m*n-m-p+1)
  lcl <- 0
  ctrl.limits <- matrix(c(lcl, ucl), ncol = 2)
  # Phase 2 prediction limits   
  ucl <- p*(m+1)*(n-1)/(m*n-m-p+1)*qf(conf, p, m*n-m-p+1)
  lcl <- 0
  pred.limits <- matrix(c(lcl, ucl), ncol = 2)
  #
  rownames(ctrl.limits) <- rownames(pred.limits) <- rep("", nrow(pred.limits))
  colnames(ctrl.limits) <- c("LCL", "UCL")
  colnames(pred.limits) <- c("LPL", "UPL")
  #
  return(list(control = ctrl.limits, prediction = pred.limits))
}

# T2 chart single observation per group



#' Statistics used in computing and drawing the Hotelling T^2 chart for
#' individual observations data
#'
#' These functions are used to compute statistics required by the \eqn{T^2}
#' chart for individual observations.
#'
#'
#' @aliases stats.T2.single limits.T2.single
#' @export stats.T2.single
#' @export limits.T2.single
#' @param data the observed data values
#' @param center a vector of values to use for center of input variables.
#' @param cov a matrix of values to use for the covariance matrix of input
#' variables.
#' @param ngroups number of groups
#' @param size sample size
#' @param nvars number of variables
#' @param conf confidence level (0 < `conf` < 1)
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
#' @references Mason, R.L. and Young, J.C. (2002) *Multivariate
#' Statistical Process Control with Industrial Applications*, SIAM.
#'
#' Montgomery, D.C. (2013) *Introduction to Statistical Quality Control*,
#' 7th ed. New York: John Wiley & Sons.
#'
#' Ryan, T. P. (2011), *Statistical Methods for Quality Improvement*, 3rd
#' ed. New York: John Wiley & Sons, Inc.
#' @keywords htest hplot multivariate
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

limits.T2.single <- function(ngroups, size = 1, nvars, conf)
{ 
  m   <- ngroups     # num. of samples
  n   <- size        # samples size
  p   <- nvars       # num. of variables
  # Phase 1 control limits
  # Tracy Mason Young (1992)
  ucl <- (m-1)^2/m*qbeta(conf, p/2, (m-p-1)/2)
  lcl <- 0
  ctrl.limits <- matrix(c(lcl, ucl), ncol = 2)
  # Phase 2 prediction limits
  ucl <- p*(m+1)*(m-1)/(m*(m-p))*qf(conf, p, m-p)
  lcl <- 0
  pred.limits <- matrix(c(lcl, ucl), ncol = 2)
  #
  rownames(ctrl.limits) <- rownames(pred.limits) <- rep("", nrow(pred.limits))
  colnames(ctrl.limits) <- c("LCL", "UCL")
  colnames(pred.limits) <- c("LPL", "UPL")
  #
  return(list(control = ctrl.limits, prediction = pred.limits))
}
