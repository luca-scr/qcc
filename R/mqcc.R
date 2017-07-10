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

mqcc <- function(data, type = c("T2", "T2.single"), center, cov,
                 limits = TRUE, pred.limits = FALSE,
                 data.name, labels, newdata, newlabels, 
                 confidence.level = (1-0.0027)^p, rules = shewhart.rules, 
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
  names(statistics) <-  labels
  names(stats$center) <- var.names
  dimnames(stats$cov) <- list(var.names, var.names)
  #
  object <- list(call = call, data.name = data.name,
                 var.names = var.names, data = data,
                 type = type,  confidence.level = confidence.level,
                 statistics = statistics,  means = stats$means,
                 center = stats$center, cov = stats$cov)
  # check for new data provided and update object
  if(!missing(newdata))
    { newdata.name <- deparse(substitute(newdata))
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
          newlabels <- seq(start+1, start+length(newstats))
          if(length(newlabels) != length(newstats))
            stop("labels must match the length of samples provided") 
        }
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
  #
  if(is.function(rules)) 
    { violations <- rules(object, run.length = 0) 
      violations$beyond.pred.limits <- 
        rules(object, run.length = 0, limits = object$pred.limits)$beyond.limits }
  else violations <- NULL
  object$violations <- violations
  #
  class(object) <- "mqcc"
  #           
  if(plot) plot(object, ...) 
  #
  return(object)
}

print.mqcc <- function(x, ...) str(x,1)

summary.mqcc <- function(object, digits =  getOption("digits"), ...)
{
  cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  data <- object$data
  m <- unique(sapply(data, nrow))     # num. of samples
  n <- unique(sapply(data, ncol))     # samples sizes
  p <- length(data)                   # num. of variables
  data.name <- object$data.name
  type <- object$type
  cat(paste(type, "chart for", data.name, "\n"))
  statistics <- object$statistics
  cat("\nSummary of group statistics:\n")
  print(summary(statistics), digits = digits, ...)
  cat("\nNumber of variables: ",p)
  cat("\nNumber of groups: ", m)
  cat("\nGroup sample size: ", n)
  cat("\n\nCenter: \n")
  print(object$center, digits = digits, ...)
  cat("\nCovariance matrix:\n")
  print(object$cov, digits = digits, ...)
  cat("|S|: ", format(det(object$cov), digits = digits), "\n")
  #
  newdata.name <- object$newdata.name
  newstats <- object$newstats
  if(!is.null(newstats)) 
    { newdata <- object$newdata
      m <- unique(sapply(newdata, nrow))     # num. of samples
      n <- unique(sapply(newdata, ncol))     # samples sizes
      p <- length(newdata)                   # num. of variables
      cat(paste("\nSummary of group statistics in ", 
                newdata.name, ":\n", sep = ""))
      print(summary(newstats), digits = digits, ...)
      cat("\nNumber of groups: ", m)
      cat("\nGroup sample size: ", n)
      cat("\n")
    }
  #
  ctrl.limits <- object$limits
  pred.limits <- object$pred.limits
  if(!is.null(ctrl.limits)) 
    { cat("\nControl limits:\n")
      .printShortMatrix(ctrl.limits, digits = digits, ...) }
  if(!is.null(pred.limits)) 
    { cat("\nPrediction limits:\n")
      .printShortMatrix(pred.limits, digits = digits, ...) }
  invisible()        
}

plot.mqcc <- function(x, add.stats = TRUE, chart.all = TRUE, 
                      label.limits = c("LCL", "UCL"),
                      label.pred.limits = c("LPL", "UPL"),
                      title, xlab, ylab, ylim, axes.las = 0,
                      digits =  getOption("digits"),
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
    { if(is.null(newstats))
         main.title <- paste(type, "Chart\nfor", data.name)
       else if(chart.all)
               main.title <- paste(type, "Chart\nfor", data.name, 
                                   "and", newdata.name)
            else main.title <- paste(type, "Chart\nfor", newdata.name) }
  else main.title <- paste(title)

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  mar <- pmax(oldpar$mar, c(4.1,4.1,3.1,2.1))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = if(add.stats) pmax(mar, c(7.6,0,0,0)) else mar)

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
       labels = if(is.null(names(statistics))) 
                   as.character(indices) else names(statistics))
  axis(2, las = axes.las)
  box()
  top.line <- par("mar")[3]-length(capture.output(cat(main.title)))
  top.line <- top.line - if(chart.all & (!is.null(newstats))) 0.1 else 0.5
  mtext(main.title, side = 3, line = top.line,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  
  lines(indices, statistics, type = "b", pch=20) 

  if(is.numeric(limits))
    { abline(h = limits, lty = 2)
      mtext(label.limits, side = 4, at = c(limits[1], limits[2]), 
           las = 1, line = 0.1, col = gray(0.3))
    }
  if(is.numeric(pred.limits))
    { abline(h = pred.limits, lty = 3)
      mtext(label.pred.limits, side = 4, 
            at = c(pred.limits[1], pred.limits[2]), 
            las = 1, line = 0.1, col = gray(0.3),
            cex = par("cex"))
    }
  #
  if(is.null(qcc.options("beyond.limits")))
     stop(".qcc.options$beyond.limits undefined. See help(qcc.options).")
  if(length(violations$beyond.limits))
    { v <- violations$beyond.limits
      if(!chart.all & !is.null(newstats))
        { v <- v - length(stats) 
          v <- v[v>0] }
      points(indices[v], statistics[v], 
             col = qcc.options("beyond.limits")$col, 
             pch = qcc.options("beyond.limits")$pch) 
    }

  if(chart.all & (!is.null(newstats)))
    { len.obj.stats <- length(object$statistics)
      len.new.stats <- length(statistics) - len.obj.stats
      abline(v = len.obj.stats + 0.5, lty = 3)
      mtext(# paste("Calibration data in", data.name),
            "Calibration data", cex = par("cex")*0.8,
            at = len.obj.stats/2, line = 0, adj = 0.5)
      mtext(# paste("New data in", object$newdata.name),  
            "New data", cex = par("cex")*0.8, 
            at = len.obj.stats + len.new.stats/2, line = 0, adj = 0.5)
  }
  #
  if(add.stats) 
    { # computes the x margins of the figure region
      plt <- par()$plt; usr <- par()$usr
      px <- diff(usr[1:2])/diff(plt[1:2])
      xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
      at.col <- xfig[1] + diff(xfig[1:2])*c(0.15, 0.45, 0.75)
      top.line <- 4.5
      # write info at bottom
      mtext(paste("Number of groups = ", length(statistics), sep = ""), 
            side = 1, line = top.line, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Sample size = ", signif(n, digits), sep = ""),
            side = 1, line = top.line+1, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("|S| = ", signif(var, digits), sep = ""),
            side = 1, line = top.line+2, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      #
      if(is.numeric(limits))
        { mtext(paste(label.limits[1], " = ", signif(limits[1], digits), sep = ""), 
                side = 1, line = top.line, adj = 0, at = at.col[2],
                font = qcc.options("font.stats"),
                cex = par("cex")*qcc.options("cex.stats"))
          mtext(paste(label.limits[2], " = ", signif(limits[2], digits), sep = ""),
                side = 1, line = top.line+1, adj = 0, at = at.col[2],
                font = qcc.options("font.stats"),
                cex = par("cex")*qcc.options("cex.stats"))
          mtext(paste("Num. beyond limits =",
                      length(unique(violations$beyond.limits))), 
                side = 1, line = top.line+2, adj = 0, at = at.col[2],
                font = qcc.options("font.stats"),
                cex = par("cex")*qcc.options("cex.stats"))
        }        
      #
      if(is.numeric(pred.limits))
        { mtext(paste(label.pred.limits[1], " = ", signif(pred.limits[1], digits), sep = ""),
                side = 1, line = top.line, adj = 0, at = at.col[3],
                font = qcc.options("font.stats"),
                cex = par("cex")*qcc.options("cex.stats"))
          mtext(paste(label.pred.limits[2], " = ", signif(pred.limits[2], digits), sep = ""),
                side = 1, line = top.line+1, adj = 0, at = at.col[3],
                font = qcc.options("font.stats"),
                cex = par("cex")*qcc.options("cex.stats"))
          mtext(paste("Num. beyond limits =",
                      length(unique(violations$beyond.pred.limits))), 
                      side = 1, line = top.line+2, adj = 0, at = at.col[3],
                font = qcc.options("font.stats"),
                cex = par("cex")*qcc.options("cex.stats"))
        }        
    }
  invisible() 
}

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
    { if(is.null(object$newstats))
           main.title <- paste("Ellipse chart\nfor", object$data.name)
      else if(chart.all)
                main.title <- paste("Ellipse chart\nfor", 
                                    object$data.name, "and",
                                    object$newdata.name)
           else main.title <- paste("Ellipse chart\nfor", 
                                    object$newdata.name) }
  else main.title <- paste(title)
  #
  grid <- cbind(seq(xlim[1], xlim[2], length = ngrid),
                seq(ylim[1], ylim[2], length = ngrid))
  x     <- grid - matrix(center, ngrid, 2, byrow=TRUE)
  cov.inv <- solve(object$cov)
  T2    <- n*apply(expand.grid(x[,1], x[,2]), 1,
                               function(x) x %*% cov.inv %*% x)
  T2    <- matrix(T2, ngrid, ngrid)
  q <- object$limits[2]
  
  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = pmax(oldpar$mar, c(4.1,4.1,3.1,2.1)))

  # plot ellipse chart
  plot(stats, type = "n", xlim = xlim, ylim = ylim, 
       ylab = if(missing(ylab)) object$var.names[2] else ylab,
       xlab = if(missing(xlab)) object$var.names[1] else xlab)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  box()
  top.line <- par("mar")[3]-length(capture.output(cat(main.title)))-0.5
  mtext(main.title, side = 3, line = top.line,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  
  contour(grid[,1], grid[,2], T2, levels = q, drawlabels = FALSE, add=TRUE)
  points(center[1], center[2], pch=3, cex=2)

  v <- object$violations$beyond.limits
  if(show.id) 
    { text(stats, labels = names(object$statistics), 
           cex = 0.8*qcc.options("cex"),
           col = { col <- rep(1, nrow(stats))
                   col[v] <- qcc.options("beyond.limits")$col
                   col })
  }
  else        
    { points(stats, ...) 
      if(length(v))
        { points(stats[v,,drop=FALSE], 
                 col = qcc.options("beyond.limits")$col, 
                 pch = qcc.options("beyond.limits")$pch)}
  }

  if(!is.null(q1) & !is.null(q2)) 
    { abline(v=q1$limits, lty=2)
      abline(h=q2$limits, lty=2) }
  #
  invisible() 
}

# T2 chart

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

stats.T2.single <- function(data, center = NULL, cov = NULL)
{ 
  data <- as.matrix(as.data.frame(data))
  m <- nrow(data)                       # num. of samples
  p <- ncol(data)                       # num. of variables
  n <- 1                                # samples sizes
  if(is.null(center))
    { center <- apply(data, 2, mean) }  # overall mean
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
