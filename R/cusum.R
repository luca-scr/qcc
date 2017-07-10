#-------------------------------------------------------------------#
#                                                                   #
#                      CUSUM CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

cusum <- function(data, sizes, center, std.dev, head.start = 0, decision.interval = 5,
                  se.shift = 1, data.name, labels, newdata, newsizes, newlabels, plot = TRUE, ...)
{
  call <- match.call()
  if (missing(data))
     stop("'data' argument is not specified")

  if(missing(data.name)) 
     data.name <- deparse(substitute(data))
  data <- data.matrix(data)

  if(missing(sizes)) 
    { sizes <- apply(data, 1, function(x) sum(!is.na(x)))  }
  else
    { if(length(sizes)==1)
         sizes <- rep(sizes, nrow(data))
      else if(length(sizes) != nrow(data))
              stop("sizes length doesn't match with data") }

  if (decision.interval <= 0)
      stop("decision.interval must be positive")

  if (head.start < 0 || head.start >= decision.interval)
      stop("head.start must be non-negative and less than decision.interval")
  
  # used for computing statistics and std.dev
  type <- if(any(sizes==1)) "xbar.one" else "xbar"
  if(ncol(data) == 1 & any(sizes > 1) & missing(std.dev))
     stop("sizes larger than 1 but data appears to be single samples. In this case you must provide also the std.dev")
  
  if(missing(labels))
    { if(is.null(rownames(data))) labels <- 1:nrow(data)
      else                        labels <- rownames(data) }

  stats <- paste("stats.", type, sep = "")
  if(!exists(stats, mode="function"))
     stop(paste("function", stats, "is not defined"))
  stats <- do.call(stats, list(data, sizes))
  statistics <- stats$statistics
  if(missing(center)) center <- stats$center

  sd <- paste("sd.", type, sep = "")
  if(!exists(sd, mode="function"))
     stop(paste("function", sd, "is not defined!"))
  if(missing(std.dev)) 
    { std.dev <- switch(type, 
                        "xbar" = { if(any(sizes > 25)) "RMSDF"
                                   else                "UWAVE-R" },
                         NULL)
      std.dev <- do.call(sd, list(data, sizes, std.dev)) }
  else 
     { if (is.character(std.dev))
          { std.dev <- do.call(sd, list(data, sizes, std.dev)) }
       else
          { if (!is.numeric(std.dev))
               stop("if provided the argument 'std.dev' must be a method available or a numerical value. See help(qcc).")  }
     }

  names(statistics) <-  rownames(data) <-  labels
  names(dimnames(data)) <- list("Group", "Samples")

  object <- list(call = call, type = "cusum", 
                 data.name = data.name, data = data, 
                 statistics = statistics, sizes = sizes, 
                 center = center, std.dev = std.dev)
  # check for new data provided and update object
  if(!missing(newdata))
    { newdata.name <- deparse(substitute(newdata))
      newdata <- data.matrix(newdata)
      if(missing(newsizes))
        { newsizes <- apply(newdata, 1, function(x) sum(!is.na(x))) }
      else
        { if(length(newsizes)==1)
            newsizes <- rep(newsizes, nrow(newdata))
          else if(length(newsizes) != nrow(newdata))
                  stop("newsizes length doesn't match with newdata") }
      stats <- paste("stats.", type, sep = "")
      if(!exists(stats, mode="function"))
         stop(paste("function", stats, "is not defined"))
      newstats <- do.call(stats, list(newdata, newsizes))$statistics
      if(missing(newlabels))
        { if(is.null(rownames(newdata)))
            { start <- length(statistics)
              newlabels <- seq(start+1, start+length(newstats)) }
          else
            { newlabels <- rownames(newdata) }
        }
      names(newstats) <- newlabels
      object$newstats <- newstats
      object$newdata  <- newdata
      object$newsizes <- newsizes
      object$newdata.name <- newdata.name
      statistics <- c(statistics, newstats)
      sizes <- c(sizes, newsizes)
    }

  n <- length(statistics)
  z <- (statistics - center)/(std.dev/sqrt(sizes))
  ldb <-  - decision.interval
  udb <- decision.interval
  #
  z.f <- z - se.shift/2
  cusum.pos <- rep(NA, n)
  cusum.pos[1] <- max(0, head.start + z.f[1])
  for (i in 2:n)
      cusum.pos[i] <- max(0, cusum.pos[i-1]+z.f[i])
  #
  z.f <- z + se.shift/2
  cusum.neg <- rep(NA, n)
  cusum.neg[1] <- max(0, head.start - z.f[1])
  for (i in 2:n)
      cusum.neg[i] <- max(0, cusum.neg[i-1]-z.f[i])
  cusum.neg <- -cusum.neg
  violations <- list(lower = which(cusum.neg < ldb),
                     upper = which(cusum.pos > udb))
  
  object$type <- "cusum"
  object$pos <- cusum.pos 
  object$neg <- cusum.neg
  object$head.start <- head.start
  object$decision.interval <- decision.interval
  object$se.shift <- se.shift
  object$violations <- violations

  class(object) <- "cusum.qcc"
  if(plot) plot(object, ...) 
  return(object)
}


print.cusum.qcc <- function(x, ...) str(x, 1)

summary.cusum.qcc <- function(object, digits =  getOption("digits"), ...)
{
  # object <- x   # Argh.  Really want to use 'object' anyway
  cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  data.name <- object$data.name
  type <- object$type
  cat(paste(type, "chart for", data.name, "\n"))
  statistics <- object$statistics
  cat("\nSummary of group statistics:\n")
  print(summary(statistics), ...)
  sizes <- object$sizes
  if(length(unique(sizes))==1)
     sizes <- sizes[1]
  if(length(sizes) == 1)
     cat("\nGroup sample size: ", format(sizes))
  else {
         cat("\nSummary of group sample sizes: ")
         tab <- table(sizes)
         print(matrix(c(as.numeric(names(tab)), tab), 
                      ncol = length(tab), byrow = TRUE, 
                      dimnames = list(c("  sizes", "  counts"),
                      character(length(tab)))),
               digits = digits, ...)
        }
  cat("\nNumber of groups: ", length(statistics))
  center <- object$center
  cat("\nCenter of group statistics: ", format(center, digits = digits))
  sd <- object$std.dev
  cat("\nStandard deviation: ", format(sd, digits = digits), "\n")

  newdata.name <- object$newdata.name
  newstats <- object$newstats
  if(!is.null(newstats)) 
    { cat(paste("\nSummary of group statistics in ", 
                newdata.name, ":\n", sep = ""))
      print(summary(newstats), digits = digits, ...)
      newsizes <- object$newsizes
      if (length(unique(newsizes)) == 1)
         newsizes <- newsizes[1]
      if(length(newsizes) == 1)
         cat("\nGroup sample size: ", format(newsizes))
      else 
         { cat("\nSummary of group sample sizes:\n")
           new.tab <- table(newsizes)
           print(matrix(c(as.numeric(names(new.tab)), new.tab),
                        ncol = length(new.tab), byrow = TRUE, 
                        dimnames = list(c("  sizes", "  counts"),
                                        character(length(new.tab)))),
                 digits = digits, ...)
         }
      cat("\nNumber of groups: ", length(newstats), "\n")
     }

  if (object$head.start > 0)
      cat("Head start (std.err.):",
          signif(object$head.start, digits = digits), "\n")
  cat("\nDecision interval (std.err.):", 
      signif(object$decision.interval, digits = digits), "\n")
  cat("Shift detection  (std. err.):", 
      signif(object$se.shift, digits = digits), "\n")

  invisible()
}


plot.cusum.qcc <- function(x, add.stats = TRUE, chart.all = TRUE, 
                           label.bounds = c("LDB", "UDB"),
                           title, xlab, ylab, ylim, axes.las = 0,
                           digits =  getOption("digits"),
                           restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "cusum.qcc")))
     stop("an object of class `cusum.qcc' is required")

  # collect info from object
  type <- object$type
  data.name <- object$data.name
  center <- object$center
  std.dev <- object$std.dev
  stats <- object$statistics
  ldb <- -object$decision.interval
  udb <- object$decision.interval
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
  cusum.pos <- object$pos
  cusum.neg <- object$neg
  beyond.bounds <- object$violations
  n.beyond.bounds <- sum(sapply(beyond.bounds, length))
  
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
  
  plot(indices, statistics, type="n",
       ylim = if(!missing(ylim)) ylim 
              else range(cusum.pos, cusum.neg, ldb, udb),
       ylab = "", xlab = if(missing(xlab)) "Group" else xlab, 
       axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  axis(1, at = indices, las = axes.las,
       labels = if(is.null(names(statistics))) 
                   as.character(indices) else names(statistics))
  axis(2, las = axes.las)
  box()
  top.line <- par("mar")[3] - length(capture.output(cat(main.title)))
  top.line <- top.line - if(chart.all & (!is.null(newstats))) 0.1 else 0.5
  mtext(main.title, side = 3, line = top.line,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))

  abline(h = 0, lwd = 2)
  abline(h = c(ldb, udb), lty = 2)
  
  lines(indices, cusum.pos[indices], type = "b", pch=20)
  lines(indices, cusum.neg[indices], type = "b", pch=20) 

  mtext(if(missing(ylab)) "Cumulative Sum" else ylab, line=3, side=2)
  lab <- "Above target"
  if (add.stats && object$head.start > 0)
      lab <- paste(lab, " (start = ", object$head.start, ")", sep = "")
  mtext(lab, srt=90, line=2, side=2, at=0+par("usr")[4]/2,
        cex = par("cex")*0.8)
  lab <- "Below target"
  if (add.stats && object$head.start > 0)
      lab <- paste(lab, " (start = ", - object$head.start, ")", sep = "")
  mtext(lab, srt=90, line=2, side=2, at=0+par("usr")[3]/2,
        cex = par("cex")*0.8)
  mtext(label.bounds, side = 4, at = c(ldb, udb), las = 1, line = 0.1, 
        col = gray(0.3), cex = par("cex"))

  if(n.beyond.bounds > 0)
    { if(is.null(qcc.options("beyond.limits")))
         stop(".qcc.options$beyond.limits undefined. See help(qcc.options).")
      v <- beyond.bounds$upper
      points(v, cusum.pos[v],
             col = qcc.options("beyond.limits")$col, 
             pch = qcc.options("beyond.limits")$pch)
      v <- beyond.bounds$lower
      points(v, cusum.neg[v], 
             col = qcc.options("beyond.limits")$col, 
             pch = qcc.options("beyond.limits")$pch)
    }
              
  if(chart.all & !is.null(newstats))
    { len.obj.stats <- length(stats)
      len.new.stats <- length(newstats)
      abline(v = len.obj.stats + 0.5, lty = 3)
      mtext(# paste("Calibration Data in", data.name),
            "Calibration data", cex = par("cex")*0.8, 
            at = len.obj.stats/2, line = 0, adj = 0.5)
      mtext(# paste("New Data in", object$newdata.name), 
            "New data", cex = par("cex")*0.8, 
            at = len.obj.stats + len.new.stats/2, line = 0, adj = 0.5)
     }

  if(add.stats) 
    { # computes the x margins of the figure region
      plt <- par()$plt; usr <- par()$usr
      px <- diff(usr[1:2])/diff(plt[1:2])
      xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
      at.col <- xfig[1] + diff(xfig[1:2])*c(0.15, 0.55)
      top.line <- 4.5
      # write info at bottom
      mtext(paste("Number of groups = ", length(statistics), sep = ""),
            side = 1, line = top.line, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      if(length(center) == 1)
         { mtext(paste("Center = ", signif(center[1], digits), sep = ""),
                 side = 1, line = top.line+1, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
         }
       else 
         { mtext("Center is variable", 
                 side = 1, line = top.line+1, adj = 0, at = at.col[1],
                 font = qcc.options("font.stats"),
                 cex = par("cex")*qcc.options("cex.stats"))
         }
      mtext(paste("StdDev = ", signif(std.dev, digits), sep = ""),
            side = 1, line = top.line+2, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Decision interval (std. err.) =", 
                  signif(object$decision.interval, digits = digits)),
            side = 1, line = top.line, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Shift detection (std. err.) =", 
                  signif(object$se.shift, digits = digits)), 
            side = 1, line = top.line+1, adj = 0, at = at.col[2], 
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("No. of points beyond boundaries =", n.beyond.bounds),
            side = 1, line = top.line+2, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      }

  invisible()
}
