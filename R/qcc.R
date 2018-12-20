#----------------------------------------------------------------------------#
#                                                                            #
#                     QUALITY CONTROL CHARTS IN R                            #
#                                                                            #
#  An R package for statistical in-line quality control.                     #
#                                                                            #
#  Written by: Luca Scrucca                                                  #
#              Department of Economics                                       #
#              University of Perugia, ITALY                                  #
#              luca.scrucca@unipg.it                                         #
#                                                                            #
#----------------------------------------------------------------------------#

#
#  Main function to create a 'qcc' object
#

qcc <- function(data, 
                type = c("xbar", "R", "S", "xbar.one", 
                         "p", "np", "c", "u", "g"),
                sizes, center, std.dev, limits, data.name, labels, 
                newdata, newsizes, newdata.name, newlabels, 
                nsigmas = 3, confidence.level, 
                rules = c(1,4), 
                plot = TRUE, ...)
{
  call <- match.call()
  
  if (missing(data))
     stop("'data' argument is not specified")
  
  if(identical(type, eval(formals(qcc)$type)))
    { type <- as.character(type)[1]
      warning("chart 'type' not specified, assuming \"", type, "\"",
              immediate. = TRUE) }
  if(!exists(paste("stats.", type, sep = ""), mode="function") |
     !exists(paste("sd.", type, sep = ""), mode="function") |
     !exists(paste("limits.", type, sep = ""), mode="function"))
    stop(paste("invalid", type, "control chart. See help(qcc) "))

  if (missing(data.name)) 
     data.name <- deparse(substitute(data))
  data <- data.matrix(data)
  if (missing(sizes)) 
     { if (any(type==c("p", "np", "u")))
          stop(paste("sample 'sizes' must be given for a", type, "Chart"))
       else
          sizes <- apply(data, 1, function(x) sum(!is.na(x)))  }
  else
     { if (length(sizes)==1)
          sizes <- rep(sizes, nrow(data))
       else if (length(sizes) != nrow(data))
                stop("sizes length doesn't match with data") }

  if (missing(labels))
     { if (is.null(rownames(data))) labels <- 1:nrow(data)
       else                         labels <- rownames(data) }

  stats <- paste("stats.", type, sep = "")
  if (!exists(stats, mode="function"))
     stop(paste("function", stats, "is not defined"))
  stats <- do.call(stats, list(data, sizes))
  statistics <- stats$statistics
  if (missing(center)) center <- stats$center

  sd <- paste("sd.", type, sep = "")
  if (!exists(sd, mode="function"))
     stop(paste("function", sd, "is not defined!"))
  missing.std.dev <- missing(std.dev)
  if (missing.std.dev)
     { std.dev <- NULL
       std.dev <- switch(type, 
                         "xbar" = { if(any(sizes > 25)) "RMSDF"
                                    else                "UWAVE-R" },
                         "xbar.one" = "MR",
                         "R" = "UWAVE-R",
                         "S" = "UWAVE-SD",
                         NULL)
       std.dev <- do.call(sd, list(data, sizes, std.dev)) }
  else 
     { if (is.character(std.dev))
          { std.dev <- do.call(sd, list(data, sizes, std.dev)) }
       else
          { if (!is.numeric(std.dev))
               stop("if provided the argument 'std.dev' must be a method available or a numerical value. See help(qcc).")  }
     }

  stopifnot(length(labels) == length(statistics))
  names(statistics) <-  rownames(data) <-  labels
  names(dimnames(data)) <- list("Group", "Samples")
  rules <- as.numeric(rules)
  
  # create object of class 'qcc'
  object <- list(call = call, type = type, rules = rules,
                 data.name = data.name, data = data, 
                 statistics = statistics, sizes = sizes, 
                 center = center, std.dev = std.dev)
  class(object) <- "qcc"

  # check for new data provided and update object
  if(!missing(newdata))
    { if (missing(newdata.name))
			{ newdata.name <- deparse(substitute(newdata))}
        newdata <- data.matrix(newdata)
        if (missing(newsizes))
           { if (any(type==c("p", "np", "u")))
                stop(paste("sample sizes must be given for a", type, "Chart"))
             else
                newsizes <- apply(newdata, 1, function(x) sum(!is.na(x))) }
        else
           { if (length(newsizes)==1)
                newsizes <- rep(newsizes, nrow(newdata))
             else if (length(newsizes) != nrow(newdata))
                      stop("newsizes length doesn't match with newdata") }
        stats <- paste("stats.", type, sep = "")
        if (!exists(stats, mode="function"))
           stop(paste("function", stats, "is not defined"))
        newstats <- do.call(stats, list(newdata, newsizes))$statistics
        if (missing(newlabels))
           { if (is.null(rownames(newdata)))
                { start <- length(statistics)
                  newlabels <- seq(start+1, start+length(newstats)) }
             else
                { newlabels <- rownames(newdata) }
        }
        stopifnot(length(newlabels) == length(newstats))
        names(newstats) <- newlabels
        object$newstats <- newstats
        object$newdata  <- newdata
        object$newsizes <- newsizes
        object$newdata.name <- newdata.name
        statistics <- c(statistics, newstats)
        sizes <- c(sizes, newsizes)
  }

  conf <- nsigmas
  if (!missing(confidence.level))
     conf <- confidence.level
  if (conf >= 1)
     { object$nsigmas <- conf }
  else
     if (conf > 0 & conf < 1)
        { object$confidence.level <- conf
          object$nsigmas <- qnorm(1 - (1 - conf)/2) }

  # get control limits
  if (missing(limits))
     { limits <- paste("limits.", type, sep = "")
       if (!exists(limits, mode="function"))
          stop(paste("function", limits, "is not defined"))
       limits <- do.call(limits, list(center = center, 
                                      std.dev = std.dev,
                                      sizes = sizes, 
                                      nsigmas = object$nsigmas,
                                      conf = object$confidence.level)) 
     }
  else 
     { if (!missing.std.dev)
          warning("'std.dev' is not used when limits is given")
       if (!is.numeric(limits))
          stop("'limits' must be a vector of length 2 or a 2-columns matrix")
       limits <- matrix(limits, ncol = 2)
       dimnames(limits) <- list(rep("",nrow(limits)), c("LCL ", "UCL"))
     }
  object$limits <- limits
  
  # get violating runs
  object$violations <- qccRules(object)

  if(plot) plot(object, ...) 
  return(object)
}

print.qcc <- function(x, digits = getOption("digits"), ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  # cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  cat(cli::rule(left = crayon::bold("Quality Control Chart"), 
                width = min(getOption("width"),50)), "\n\n")

  data.name <- object$data.name
  type <- object$type
  statistics <- object$statistics
  # cat(paste(type, "chart for", data.name, "\n"))
  # cat("\nSummary of group statistics:\n")
  # print(summary(statistics), digits = digits, ...)
  
  cat("Chart type                 =", type, "\n")
  cat("Data (phase I)             =", data.name, "\n")
  cat("Number of groups           =", length(statistics), "\n")

  sizes <- object$sizes
  if(length(unique(sizes))==1)
     sizes <- sizes[1]
  if(length(sizes) == 1)
  {
    cat("Group sample size          =", signif(sizes), "\n")
  } else 
  {
    cat("Group sample sizes         =")
    tab <- table(sizes)
    print(matrix(c(as.numeric(names(tab)), tab), 
                 ncol = length(tab), byrow = TRUE, 
                 dimnames = list(c("  sizes", "  counts"),
                                 character(length(tab)))), 
          digits = digits, ...)
  }
  
  center <- object$center
  if(length(center) == 1)
  { 
    cat("Center of group statistics =", signif(center, digits = digits), "\n")   } else
  { 
    out <- paste(signif(center, digits = digits))
    out <- out[which(cumsum(nchar(out)+1) < getOption("width")-40)]      
    out <- paste0(paste(out, collapse = " "), " ...")
    cat("Center of group statistics =", out, "\n", sep = "")
  }
  
  sd <- object$std.dev
  if(length(sd) == 1)
  { 
    cat("Standard deviation         =", signif(sd, digits = digits), "\n") 
  } else
  { 
    out <- paste(signif(sd, digits = digits))
    out <- out[which(cumsum(nchar(out)+1) < getOption("width")-40)]
    out <- paste0(paste(out, collapse = " "), " ...")
    cat("Standard deviation         =", out, "\n", sep = "")
  }

  newdata.name <- object$newdata.name
  newstats <- object$newstats
  if (!is.null(newstats)) 
  { 
    # cat(cli::rule(line = 1, width = 50), "\n")
    # cat(paste("\nSummary of group statistics in ", 
    #           newdata.name, ":\n", sep = ""))
    # print(summary(newstats), digits = digits, ...)
    cat("\nNew data (phase II)        =", newdata.name, "\n")
    cat("Number of groups           =", length(newstats), "\n")
    newsizes <- object$newsizes
    if (length(unique(newsizes)) == 1)
      newsizes <- newsizes[1]
    if (length(newsizes) == 1)
    {
      cat("Group sample size          =", signif(newsizes), "\n")
    } else 
    { cat("Group sample sizes         =")
      new.tab <- table(newsizes)
      print(matrix(c(as.numeric(names(new.tab)), new.tab),
                   ncol = length(new.tab), byrow = TRUE, 
                   dimnames = list(c("  sizes", "  counts"),
                                   character(length(new.tab)))), 
            digits = digits, ...)
    }
  }
  
  # cat(cli::rule(line = 1, width = 50), "\n")
  limits <- object$limits
  if(!is.null(limits)) 
  { 
    # cat("Control limits:\n")
    cat("\nControl limits at nsigmas  =", object$nsigmas, "\n")    
    # names(dimnames(limits)) <- c("Control limits             =", "")
    .printShortMatrix(limits, digits = digits, ...) 
  }

  invisible()
}

summary.qcc <- function(object, ...) print.qcc(object, ...)


plot.qcc <- function(x, 
                     add.stats = qcc.options("add.stats"), 
                     chart.all = qcc.options("chart.all"), 
                     fill = qcc.options("fill"),
                     label.center = "CL",
                     label.limits = c("LCL ", "UCL"), 
                     title, xlab, ylab, ylim, axes.las = 0,
                     digits = getOption("digits"),
                     restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "qcc")))
    stop("an object of class `qcc' is required")
  
  # collect info from object
  type <- object$type
  std.dev <- object$std.dev
  data.name <- object$data.name
  center <- object$center
  stats <- object$statistics
  limits <- object$limits 
  lcl <- limits[,1]
  ucl <- limits[,2]
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
  if(chart.all) 
    { statistics <- c(stats, newstats)
      indices <- 1:length(statistics) 
    } else 
    { if(is.null(newstats))
        { statistics <- stats
          indices <- 1:length(statistics) 
        } else 
        { statistics <- newstats 
          indices <- seq(length(stats)+1, length(stats)+length(newstats)) 
        }
    }
  
  if(missing(title))
  { 
    if(is.null(newstats))
      main.title <- paste(type, "chart for", data.name)
    else if(chart.all)
           main.title <- paste(type, "chart for", data.name, 
                                "and", newdata.name)
         else 
           main.title <- paste(type, "chart for", newdata.name) 
  }
  else main.title <- paste(title)
  cex.labels <- par("cex")*qcc.options("cex")
  cex.stats <- par("cex")*qcc.options("cex.stats")
  
  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      # mgp = c(2.1, 0.8, 0),
      mar = pmax(par("mar"), c(4.1,4.1,1.1,2.1), na.rm=TRUE),
      oma = if(add.stats) c(3.5*cex.stats, 0, 1.5*cex.labels, 0) 
            else          c(0, 0, 1.5*cex.labels, 0))
  
  # plot Shewhart chart
  plot(indices, statistics, type="n",
       ylim = if(!missing(ylim)) ylim 
              else range(statistics, limits, center),
       ylab = if(missing(ylab)) "Group summary statistics" else ylab,
       xlab = if(missing(xlab)) "Group" else xlab, 
       axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  axis(1, at = indices, las = axes.las,
       labels = if(is.null(names(statistics))) 
                as.character(indices) else names(statistics),
       cex.axis = par("cex.axis")*0.9)
  axis(2, las = axes.las, cex.axis = par("cex.axis")*0.9)
  box()
  mtext(main.title, side = 3, outer = TRUE, 
        line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = par("cex")*qcc.options("cex"), 
        col  = par("col.main"))

  # draw control limits
  if(any(object$rules == 1))
  { if(fill)
    { # fill the in-control area
      if(length(lcl) == 1) 
      { 
        polygon(par("usr")[c(1,2,2,1)], c(lcl, lcl, ucl, ucl), border = FALSE, 
                col = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.1))
      } else
      {
        x1 <- x2 <- c(indices, indices[length(indices)]+1)-0.5
        y1 <- c(lcl[indices], lcl[indices[length(indices)]])
        y2 <- c(ucl[indices], ucl[indices[length(indices)]])
        xp1 <- rep(x1, each=2)[-1]
        xp2 <- rep(x2, each=2)[-1]
        yp1 <- rep(y1, each=2)[-2*length(y1)]
        yp2 <- rep(y2, each=2)[-2*length(y2)]
        # plot(indices, statistics, ylim = range(statistics, limits, center))
        # lines(x1, y1, type = "s", lty = 2)
        # lines(x2, y2, type = "s", lty = 2)
        polygon(c(xp1,rev(xp2)), c(yp1,rev(yp2)), border = FALSE, 
                col = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.1))
        # x <- c(1:5, 1:5)
        # y <- c(2,1,1,3,2,5,6,6,4,5)
        # x1 = c(x[1:5],x[5]+1)-0.5
        # y1 = c(y[1:5],y[5])
        # x2 = c(x[6:10],x[10]+1)-0.5
        # y2 = c(y[6:10],y[10])
        # xp1 <- rep(x1, each=2)[-1]
        # xp2 <- rep(x2, each=2)[-1]
        # yp1 <- rep(y1, each=2)[-2*length(y1)]
        # yp2 <- rep(y2, each=2)[-2*length(y2)]
        # plot(x, y, type = "p", ylim = range(y), xlim = range(x), pch = 20)
        # lines(x1, y1, type = "s", lty = 2)
        # lines(x2, y2, type = "s", lty = 2)
        # polygon(c(xp1,rev(xp2)), c(yp1,rev(yp2)), border = FALSE, col = qcc.options("fill"))
      }
  } else
    {
      if(length(lcl) == 1) 
      { 
        abline(h = c(lcl,ucl),
               lty = qcc.options("zones")$lty[1], 
               col = qcc.options("zones")$col[1])
      } else 
      { 
        x1 <- x2 <- c(indices, indices[length(indices)]+1)-0.5
        y1 <- c(lcl[indices], lcl[indices[length(indices)]])
        y2 <- c(ucl[indices], ucl[indices[length(indices)]])
        lines(x1, y1, type="s", 
              lty = qcc.options("zones")$lty[1],
              col = qcc.options("zones")$col[1])
        lines(x2, y2, type="s", 
              lty = qcc.options("zones")$lty[1],
              col = qcc.options("zones")$col[1]) 
      }
  }
  mtext(label.limits, side = 4, at = c(rev(lcl)[1], rev(ucl)[1]), 
        las = 1, line = 0.1, col = gray(0.3), 
        cex = par("cex")*qcc.options("cex.stats"))
  mtext(label.center, side = 4, at = rev(center)[1], 
        las = 1, line = 0.1, col = gray(0.3), 
        cex = par("cex")*qcc.options("cex.stats"))
  }
  
  # draw 2-sigma warning limits
  if(any(object$rules == 2))
  { 
    limits.2sigma <- do.call(paste("limits.", object$type, sep = ""), 
                             list(center = object$center, 
                                  std.dev = object$std.dev,
                                  sizes = c(object$sizes, object$newsizes), 
                                  nsigmas = object$nsigmas*2/3))
    if(fill)
    { 
      if(nrow(limits.2sigma) == 1) 
      { 
        polygon(par("usr")[c(1,2,2,1)], limits.2sigma[1,c(1,1,2,2)],
                border = FALSE, 
                col = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.1))
      } else
      {
        x1 <- x2 <- c(indices, indices[length(indices)]+1)-0.5
        y1 <- c(limits.2sigma[indices,1], 
                limits.2sigma[indices[length(indices)],1])
        y2 <- c(limits.2sigma[indices,2], 
                limits.2sigma[indices[length(indices)],2])
        xp1 <- rep(x1, each=2)[-1]
        xp2 <- rep(x2, each=2)[-1]
        yp1 <- rep(y1, each=2)[-2*length(y1)]
        yp2 <- rep(y2, each=2)[-2*length(y2)]
        polygon(c(xp1,rev(xp2)), c(yp1,rev(yp2)), 
                border = FALSE, 
                col = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.1))
      }
    } else
    {
      if(nrow(limits.2sigma) == 1) 
      { 
        abline(h = limits.2sigma, 
               lty = qcc.options("zones")$lty[2], 
               col = qcc.options("zones")$col[2])
      } else
      { 
        x1 <- x2 <- c(indices, indices[length(indices)]+1)-0.5
        y1 <- c(limits.2sigma[indices,1], 
                limits.2sigma[indices[length(indices)],1])
        y2 <- c(limits.2sigma[indices,2], 
                limits.2sigma[indices[length(indices)],2])
        lines(x1, y1, type="s", 
              lty = qcc.options("zones")$lty[2],
              col = qcc.options("zones")$col[2])
        lines(x2, y2, type="s", 
              lty = qcc.options("zones")$lty[2],
              col = qcc.options("zones")$col[2])  
      }
    }
  }
  
  # draw 1-sigma warning limits
  if(any(object$rules == 3))
  { 
    limits.1sigma <- do.call(paste("limits.", object$type, sep = ""), 
                             list(center = object$center, 
                                  std.dev = object$std.dev,
                                  sizes = c(object$sizes, object$newsizes), 
                                  nsigmas = object$nsigmas*1/3))
    if(fill)
    { 
      if(nrow(limits.1sigma) == 1) 
      { 
        polygon(par("usr")[c(1,2,2,1)], limits.1sigma[1,c(1,1,2,2)],
                border = FALSE, 
                col = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.1))
      } else
      {
        x1 <- x2 <- c(indices, indices[length(indices)]+1)-0.5
        y1 <- c(limits.1sigma[indices,1], 
                limits.1sigma[indices[length(indices)],1])
        y2 <- c(limits.1sigma[indices,2], 
                limits.1sigma[indices[length(indices)],2])
        xp1 <- rep(x1, each=2)[-1]
        xp2 <- rep(x2, each=2)[-1]
        yp1 <- rep(y1, each=2)[-2*length(y1)]
        yp2 <- rep(y2, each=2)[-2*length(y2)]
        polygon(c(xp1,rev(xp2)), c(yp1,rev(yp2)), 
                border = FALSE, 
                col = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.1))
      }
    } else
    {
      if(nrow(limits.1sigma) == 1) 
      { 
        abline(h = limits.1sigma, 
               lty = qcc.options("zones")$lty[3], 
               col = qcc.options("zones")$col[3])
      } else
      { 
        x1 <- x2 <- c(indices, indices[length(indices)]+1)-0.5
        y1 <- c(limits.1sigma[indices,1], 
                limits.1sigma[indices[length(indices)],1])
        y2 <- c(limits.1sigma[indices,2], 
                limits.1sigma[indices[length(indices)],2])
        lines(x1, y1, type="s", 
              lty = qcc.options("zones")$lty[3],
              col = qcc.options("zones")$col[3])
        lines(x2, y2, type="s", 
              lty = qcc.options("zones")$lty[3],
              col = qcc.options("zones")$col[3])  
      }
    }
  }
  
  # if(any(object$rules == 3) & !fill)
  # {
  #   limits.1sigma <- do.call(paste("limits.", object$type, sep = ""), 
  #                            list(center = object$center, 
  #                                 std.dev = object$std.dev,
  #                                 sizes = c(object$sizes, object$newsizes), 
  #                                 nsigmas = object$nsigmas*1/3))
  #   if(nrow(limits.1sigma) == 1) 
  #   { 
  #     abline(h = limits.1sigma, 
  #            lty = qcc.options("zones")$lty[3], 
  #            col = qcc.options("zones")$col[3])
  #   } else
  #   { 
  #     x1 <- x2 <- c(indices, indices[length(indices)]+1)-0.5
  #     y1 <- c(limits.1sigma[indices,1], limits.1sigma[indices[length(indices)],1])
  #     y2 <- c(limits.1sigma[indices,2], limits.1sigma[indices[length(indices)],2])
  #     lines(x1, y1, type="s", 
  #           lty = qcc.options("zones")$lty[3],
  #           col = qcc.options("zones")$col[3])
  #     lines(x2, y2, type="s", 
  #           lty = qcc.options("zones")$lty[3],
  #           col = qcc.options("zones")$col[3])  
  #   }
  # }
  # 
  # draw center line
  if(length(center) == 1)
    abline(h = center, col = qcc.options("zones")$col[1])
  else 
    lines(indices, center[indices], type="s", 
          col = qcc.options("zones")$col[1])
  
  # draw lines & points
  lines(indices, statistics, type = "b", pch=NA)
  col <- rep(palette()[1], length(indices))
  pch <- rep(20, length(indices))
  
  if(!is.null(violations))
  { 
    for(j in 1:4)
       {
         i <- (indices %in% which(violations==j))
         col[i] <- qcc.options("rules")$col[j]
         pch[i] <- qcc.options("rules")$pch[j]
       }
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
    # seq(par("plt")[1], par("plt")[2], length.out = 4)[-4]
    # write info at bottom
    mtext(paste("Number of groups = ", length(statistics), sep = ""), 
          side = 1, outer = TRUE, line = 0, adj = 0, at = at[1],
          font = qcc.options("font.stats"),
          cex = par("cex")*qcc.options("cex.stats"))
    center <- object$center
    if(length(center) == 1)
    { mtext(paste("Center = ", signif(center[1], digits), sep = ""),
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    } else 
    { mtext("Center is variable",
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    
    if(length(std.dev) == 1)
    { mtext(paste("StdDev = ", signif(std.dev, digits), sep = ""),
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    } else
    { mtext("StdDev is variable",
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    
    if(length(unique(lcl)) == 1)
    { mtext(paste("LCL = ", signif(lcl[1], digits), sep = ""), 
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    else 
    { mtext("LCL is variable", 
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    
    if(length(unique(ucl)) == 1)
    { mtext(paste("UCL = ", signif(ucl[1], digits), sep = ""),
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats")) 
    }
    else 
    { mtext("UCL is variable", 
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
    
    if(!is.null(violations))
    { mtext(paste("Number beyond limits =",
                  sum(violations==1, na.rm=TRUE)),
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Number violating runs =",
                  sum(violations > 1, na.rm=TRUE)),
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }
  }
  
  invisible()
}

#
#  Functions used to compute Shewhart charts statistics
#

qcc.c4 <- function(n)
{ sqrt(2/(n - 1)) * exp(lgamma(n/2) - lgamma((n - 1)/2)) }

# xbar

stats.xbar <- function(data, sizes)
{
  data <- as.matrix(data)
  if(missing(sizes))
    sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  statistics <- apply(data, 1, mean, na.rm=TRUE)
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.xbar <- function(data, sizes, std.dev = c("UWAVE-R", "UWAVE-SD", "MVLUE-R", "MVLUE-SD", "RMSDF"), ...)
{
  data <- as.matrix(data)
  if(missing(sizes))
    sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  if(any(sizes == 1))
    stop("group sizes must be larger than one")
  c4 <- qcc.c4
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
                            sd <- sum(S/c4(sizes))/length(sizes) 
                          },
             "MVLUE-R"  = { R <- apply(data, 1, function(x) 
                            diff(range(x, na.rm = TRUE)))
                            d2 <- qcc.options("exp.R.unscaled")[sizes]
                            d3 <- qcc.options("se.R.unscaled")[sizes]
                            w  <- (d2/d3)^2
                            sd <- sum(R/d2*w)/sum(w) 
                          }, 
             "MVLUE-SD" = { S <- apply(data, 1, sd, na.rm = TRUE)
                            w  <- c4(sizes)^2/(1-c4(sizes)^2)
                            sd <- sum(S/c4(sizes)*w)/sum(w) 
                          },
             "RMSDF" =    { S <- apply(data, 1, sd, na.rm = TRUE)
                            w  <- sizes-1
                            sd <- sqrt(sum(S^2*w)/sum(w))/c4(sum(w)+1) 
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
  if(is.null(conf))
     { lcl <- center - nsigmas * se.stats
       ucl <- center + nsigmas * se.stats
     }
  else 
     { if (conf > 0 & conf < 1) 
          { nsigmas <- qnorm(1 - (1 - conf)/2)
            lcl <- center - nsigmas * se.stats
            ucl <- center + nsigmas * se.stats
          }
       else stop("invalid 'conf' argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}


# S chart

stats.S <- function(data, sizes)
{
  data <- as.matrix(data)
  if (missing(sizes))
     sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  if(ncol(data)==1) 
    { statistics <- as.vector(data) }
  else 
    { statistics <- sqrt(apply(data, 1, var, na.rm=TRUE)) }
  if (length(sizes == 1))
     sizes <- rep(sizes, length(statistics))
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.S <- function(data, sizes, std.dev = c("UWAVE-SD", "MVLUE-SD", "RMSDF"), ...)
{
  if (!is.numeric(std.dev))
     std.dev <- match.arg(std.dev)
  sd.xbar(data, sizes, std.dev)
}

limits.S <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  if(length(unique(sizes))==1) sizes <- sizes[1]
  c4 <- qcc.c4
  se.stats <- std.dev * sqrt(1 - c4(sizes)^2)
  if (is.null(conf)) 
     { lcl <- pmax(0, center - nsigmas * se.stats)
       ucl <- center + nsigmas * se.stats
     }
  else 
     { if (conf > 0 & conf < 1) 
          { ucl <- std.dev * sqrt(qchisq(1 - (1 - conf)/2, sizes - 1)/
                                  (sizes - 1))
            lcl <- std.dev * sqrt(qchisq((1 - conf)/2, sizes - 1)/
                                  (sizes - 1))
          }
          else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  limits
}

# R Chart 

stats.R <- function(data, sizes)
{
  data <- as.matrix(data)
  if (missing(sizes))
     sizes <- apply(data, 1, function(x) sum(!is.na(x)))
  if(ncol(data)==1) 
    { statistics <- as.vector(data) }
  else 
    { statistics <- apply(data, 1, function(x) diff(range(x, na.rm=TRUE))) }
  if (length(sizes == 1))
     sizes <- rep(sizes, length(statistics))
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.R <- function(data, sizes, std.dev = c("UWAVE-R", "MVLUE-R"), ...)
{
  if (!is.numeric(std.dev))
     std.dev <- match.arg(std.dev)
  sd.xbar(data, sizes, std.dev)
}

limits.R <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  if (length(unique(sizes))==1) sizes <- sizes[1]
  se.R.unscaled <- qcc.options("se.R.unscaled")
  Rtab <- length(se.R.unscaled)
  if (is.null(conf)) 
     { if (any(sizes > Rtab))
          stop(paste("group size must be less than", 
                      Rtab + 1, "when giving nsigmas"))
       se.R <- se.R.unscaled[sizes] * std.dev
       lcl <- pmax(0, center - nsigmas * se.R)
       ucl <- center + nsigmas * se.R
     }
  else 
     { if (conf > 0 && conf < 1) 
          { ucl <- qtukey(1 - (1 - conf)/2, sizes, 1e100) * std.dev
            lcl <- qtukey((1 - conf)/2, sizes, 1e100) * std.dev
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# xbar Chart for one-at-time data

stats.xbar.one <- function(data, sizes)
{
  statistics <- as.vector(data)
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

sd.xbar.one <- function(data, sizes, std.dev = c("MR", "SD"), r = 2, ...)
{
  data <- as.vector(data)
  n <- length(data)
  if(!is.numeric(std.dev)) 
     std.dev <- match.arg(std.dev)
  c4 <- qcc.c4
  if(is.numeric(std.dev)) 
    { sd <- std.dev }
  else
    { switch(std.dev, 
             "MR" = { d2 <- qcc.options("exp.R.unscaled")
                      if(is.null(d2))
                         stop(".qcc.options$exp.R.unscaled is null")
                      d <- 0
                      for(j in r:n)
                          d <- d+abs(diff(range(data[c(j:(j-r+1))])))
                      sd <- (d/(n-r+1))/d2[r] },
             "SD" = { sd <- sd(data)/c4(n) },
             sd <- NULL)
    }
  return(sd)
}


limits.xbar.one <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  se.stats <- std.dev
  if (is.null(conf)) 
     { lcl <- center - nsigmas * se.stats
       ucl <- center + nsigmas * se.stats
     }
  else 
     { if (conf > 0 & conf < 1) 
          { nsigmas <- qnorm(1 - (1 - conf)/2)
            lcl <- center - nsigmas * se.stats
            ucl <- center + nsigmas * se.stats
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}


# p Chart

stats.p <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  list(statistics = data/sizes, center = pbar)
}

sd.p <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  std.dev <- sqrt(pbar * (1 - pbar))
  return(std.dev)
}

limits.p <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{ 
  limits.np(center * sizes, std.dev, sizes, nsigmas, conf) / sizes
}

# np Chart

stats.np <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  center <- sizes * pbar
  if (length(unique(center)) == 1)
     center <- center[1]
  list(statistics = data, center = center)
}

sd.np <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  std.dev <- sqrt(sizes * pbar * (1 - pbar))
  if (length(unique(std.dev)) == 1)
     std.dev <- std.dev[1]
  return(std.dev)
}

limits.np <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{ 
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) == 1) sizes <- sizes[1]
  pbar <- mean(center / sizes)
  if (is.null(conf))
     { tol <- nsigmas * sqrt(pbar * (1 - pbar) * sizes)
       lcl <- pmax(center - tol, 0)
       ucl <- pmin(center + tol, sizes)
     }
  else
     { if (conf > 0 & conf < 1)
          { lcl <- qbinom((1 - conf)/2, sizes, pbar)
            ucl <- qbinom((1 - conf)/2, sizes, pbar, lower.tail = FALSE)
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# c Chart

stats.c <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  if (length(unique(sizes)) != 1)
     stop("all sizes must be be equal for a c chart")
  statistics <- data
  center <- mean(statistics)
  list(statistics = statistics, center = center)
}

sd.c <- function(data, sizes, ...)
{
  data <- as.vector(data)
  std.dev <- sqrt(mean(data))
  return(std.dev)
}

limits.c <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  if (is.null(conf))
     { lcl <- center - nsigmas * sqrt(center)
       lcl[lcl < 0] <- 0
       ucl <- center + nsigmas * sqrt(center)
     }
  else 
     { if (conf > 0 & conf < 1) 
          { ucl <- qpois(1 - (1 - conf)/2, center)
            lcl <- qpois((1 - conf)/2, center)
          }
       else stop("invalid conf argument. See help.")
     }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# u Chart

stats.u <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  statistics <- data/sizes
  center <- sum(sizes * statistics)/sum(sizes)
  list(statistics = statistics, center = center)
}

sd.u <- function(data, sizes, ...)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  std.dev <- sqrt(sum(data)/sum(sizes))
  return(std.dev)
}

limits.u <- function(center, std.dev, sizes, nsigmas = NULL, conf = NULL)
{
  if(is.null(nsigmas) & is.null(conf))
    stop("Argument 'nsigmas' or 'conf' must be provided. See help.")
  sizes <- as.vector(sizes)
  if (length(unique(sizes))==1) sizes <- sizes[1]
  limits.c(center * sizes, std.dev, sizes, nsigmas, conf) / sizes
}

