#-------------------------------------------------------------------#
#                                                                   #
#                       EWMA CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

ewmaSmooth <- function(x, y, lambda = 0.20, start, ...)
{
#
# Exponential-Weighted Moving Average 
# 
# Return smooth values based on 
# 
# z_t = lambda*y_t + (1-lambda)*z_t-1      
# 
# where 0<= lambda <=1 is the parameter which controls the weights applied 
# to the data, and start is the starting value.
# Returns a list with elements:
# x = ordered x-values
# y = smoothed fitted values of y
# 
  if (length(y)!=length(x))
     stop("x and y must have the same length!")
  if (abs(lambda)>1)
     stop("lambda parameter must be between 0 and 1")
  ord <- order(x) 
  x <- x[ord]
  y <- y[ord]
  n <- length(y)
  if (missing(start)) start <- y[1]
  z <- c(start, y)
  for (i in 2:(n + 1))
    z[i] <- lambda * z[i] + (1 - lambda) * z[i - 1]
  list(x=x, y=z[-1], lambda=lambda, start=start)
}


ewma <- function(data, 
                 sizes, center, std.dev, 
                 lambda = 0.2, nsigmas = 3, 
                 newdata, newsizes, ...)
{

  call <- match.call()
  if (missing(data))
     stop("'data' argument is not specified")

  data.name <- deparse(substitute(data))
  data <- data.matrix(data)

  if(missing(sizes)) 
    { sizes <- apply(data, 1, function(x) sum(!is.na(x)))  }
  else
    { if(length(sizes)==1)
         sizes <- rep(sizes, nrow(data))
      else if(length(sizes) != nrow(data))
              stop("sizes length doesn't match with data") }
  # used for computing statistics and std.dev
  type <- if(any(sizes==1)) "xbar.one" else "xbar"

  labels <- if(is.null(rownames(data))) 1:nrow(data) else rownames(data)

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

  stopifnot(length(labels) == length(statistics))
  names(statistics) <- rownames(data) <-  labels
  names(dimnames(data)) <- list("Group", "Samples")
  object <- list(call = call, type = "ewma", 
                 data.name = data.name, data = data, 
                 statistics = statistics, sizes = sizes, 
                 center = center, std.dev = std.dev)

  # check for new data provided and update object
  if(!missing(newdata))
  { 
    newdata.name <- deparse(substitute(newdata))
    newdata <- data.matrix(newdata)
    if(missing(newsizes))
    { 
      newsizes <- apply(newdata, 1, function(x) sum(!is.na(x))) 
    } else
    { 
      if(length(newsizes)==1)
        newsizes <- rep(newsizes, nrow(newdata))
      else 
        if(length(newsizes) != nrow(newdata))
           stop("newsizes length doesn't match with newdata") 
    }
    stats <- paste("stats.", type, sep = "")
    if(!exists(stats, mode="function"))
      stop(paste("function", stats, "is not defined"))
    newstats <- do.call(stats, list(newdata, newsizes))$statistics
    if(is.null(rownames(newdata)))
    { 
      start <- length(statistics)
      newlabels <- seq(start+1, start+length(newstats)) 
    } else
    { 
      newlabels <- rownames(newdata)
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

  n <- length(statistics)
  indices <- 1:length(statistics)
  ewma <- ewmaSmooth(indices, statistics, lambda=lambda, start=center)
  sigma2 <- std.dev^2/sizes * 
            ((lambda/(2-lambda))*(1-(1-lambda)^(2*(1:n))))
  ucl <- center + nsigmas*sqrt(sigma2)
  lcl <- center - nsigmas*sqrt(sigma2)

  object$x <- ewma$x
  y <- as.vector(ewma$y)
  names(y) <- c(names(object$statistics), names(object$newstats))
  object$ewma <- y
  object$sigma <- sqrt(sigma2)
  object$lambda <- lambda
  object$nsigmas <- nsigmas
  limits <- cbind(lcl,ucl)
  colnames(limits) <- c("LCL", "UCL")
  object$limits <- limits
  object$violations <- ifelse(y < lcl | y > ucl, 1, NA)

  class(object) <- "ewma.qcc"
  return(object)
}


print.ewma.qcc <- function(x, digits =  getOption("digits"), ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  # cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  cat(cli::rule(left = crayon::bold("EWMA Chart"), 
                width = min(getOption("width"),50)), "\n\n")
    
  data.name <- object$data.name
  # type <- object$type
  statistics <- object$statistics
  # cat("\nSummary of group statistics:\n")
  # print(summary(statistics), digits = digits, ...)

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
  cat("Center of group statistics =", signif(center, digits = digits), "\n")

  sd <- object$std.dev
  cat("Standard deviation         =", signif(sd, digits = digits), "\n")

  newdata.name <- object$newdata.name
  newstats <- object$newstats
  if(!is.null(newstats)) 
  { 
    # cat(paste("\nSummary of group statistics in ", 
    #           newdata.name, ":\n", sep = ""))
    # print(summary(newstats), digits = digits, ...)
    cat("\nNew data (phase II)        =", newdata.name, "\n")
    cat("Number of groups           =", length(newstats), "\n")
    newsizes <- object$newsizes
    if (length(unique(newsizes)) == 1)
      newsizes <- newsizes[1]
    if(length(newsizes) == 1)
    {
      cat("Group sample size          =", signif(newsizes), "\n")
    } else 
    { 
      cat("Group sample sizes         =")
      new.tab <- table(newsizes)
      print(matrix(c(as.numeric(names(new.tab)), new.tab),
                   ncol = length(new.tab), byrow = TRUE, 
                   dimnames = list(c("  sizes", "  counts"),
                                   character(length(new.tab)))),
            digits = digits, ...)
    }
  }

  cat("\nSmoothing parameter        =", 
      signif(object$lambda, digits = digits), "\n")

  limits <- object$limits
  if(!is.null(limits)) 
  { 
    cat("Control limits at nsigmas  =", object$nsigmas, "\n")    
    # names(dimnames(limits)) <- c("Control limits             =", "")
    .printShortMatrix(limits, digits = digits, ...) 
  }

  invisible()
}

summary.ewma.qcc <- function(object, ...) print.ewma.qcc(object, ...)


plot.ewma.qcc <- function(x, xtime = NULL,
                          add.stats = qcc.options("add.stats"), 
                          chart.all = qcc.options("chart.all"), 
                          fill = qcc.options("fill"),
                          label.center = "CL",
                          label.limits = c("LCL", "UCL"), 
                          title, xlab, ylab, xlim, ylim,
                          digits = getOption("digits"), 
                          ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "ewma.qcc")))
     stop("an object of class `ewma.qcc' is required")

  # collect info from object
  type <- object$type
  data.name <- object$data.name
  center <- object$center
  std.dev <- object$std.dev
  stats <- object$statistics
  ewma <- object$ewma
  limits <- object$limits
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
  statistics <- c(stats, newstats)
  groups <- if(is.null(xtime)) 1:length(statistics) else xtime
  stopifnot(length(groups) == length(statistics))

  if(missing(title))
  { 
    if(is.null(newstats))
       title <- paste(type, "chart for", data.name)
    else if(chart.all)
           title <- paste(type, "chart for", data.name, "and", newdata.name)
         else 
           title <- paste(type, "Chart for", newdata.name) 
  }

  df <- data.frame(group = groups, 
                   stat = statistics,
                   ewma = ewma,
                   limits = limits,
                   violations = factor(ifelse(is.na(violations), 0, violations),
                                             levels = 0:1),
                   check.names = FALSE)
  if(!chart.all & (!is.null(newstats)))
    df <- df[df$group > length(object$statistics),]
  
  if(missing(ylim))
    ylim <- range(df[,c("stat", "limits.LCL", "limits.UCL")], na.rm = TRUE)
  if(missing(xlim))
    xlim <- range(df$group, na.rm = TRUE)
  
  plot <- 
    ggplot(data = df, aes_string(x = "group", y = "ewma")) +
    geom_line() +
    geom_point(aes_string(colour = "violations", 
                          shape = "violations"), 
               size = 2) +
    scale_colour_manual(values = c("black", qcc.options("rules")$col)) +
    scale_shape_manual(values = c(20, qcc.options("rules")$pch)) +
    geom_point(aes_string(y = "stat"), pch = 3) +
    labs(title = title, subtitle = "",
         x = if(missing(xlab)) "Group" else xlab,
         y = if(missing(ylab)) "Group Summary Statistics" else ylab) +
    coord_cartesian(xlim = xlim+c(-0.5,0.5), 
                    ylim = extendrange(ylim),
                    expand = FALSE, clip = "off") +
    theme_light() + 
    theme(plot.background = element_rect(fill = qcc.options("bg.margin"),
                                         color = qcc.options("bg.margin")),
          panel.background = element_rect(fill = qcc.options("bg.figure")),
          plot.title = element_text(face = "bold", size = 11),
          legend.position = "none",
          plot.margin = margin(5, 30, 5, 5),
          axis.text.y = element_text(angle = 90, 
                                     margin = margin(l = 5, r = 5),
                                     hjust = 0.5, vjust = 0.5))
  
  plot <- plot + 
  {
    if(is.numeric(df$group))
      scale_x_continuous(breaks = pretty(df$group, n = 7))
    else
      scale_x_date(breaks = pretty(df$group, n = 7))
  }
   
  # draw control limits
  if(all(is.finite(limits)))
  { 
    x1 <- x2 <- c(df$group, df$group[length(df$group)]+1)-0.5
    y1 <- c(df$limits.LCL, df$limits.LCL[length(df$group)])
    y2 <- c(df$limits.UCL, df$limits.UCL[length(df$group)])
    xp1 <- rep(x1, each=2)[-1]
    xp2 <- rep(x2, each=2)[-1]
    yp1 <- rep(y1, each=2)[-2*length(y1)]
    yp2 <- rep(y2, each=2)[-2*length(y2)]
    if(fill)
    { # fill the in-control area
      plot <- plot + 
        geom_polygon(data = data.frame(x = c(xp1,rev(xp2)), 
                                       y = c(yp1,rev(yp2))),
                     aes_string(x = "x", y = "y"), 
                     fill = adjustcolor(qcc.options("zones")$fill, alpha.f=0.2),
                     col = NA)
    } else
    {
      plot <- plot + 
        geom_step(data = data.frame(x = x1, y = y1),
                  aes_string(x = "x", y = "y"), 
                  lty = qcc.options("zones")$lty[1],
                  col = qcc.options("zones")$col[1])
      plot <- plot + 
        geom_step(data = data.frame(x = x2, y = y2),
                  aes_string(x = "x", y = "y"), 
                  lty = qcc.options("zones")$lty[1],
                  col = qcc.options("zones")$col[1])
    }

    plot <- plot + 
      geom_text(data = data.frame(y = c(rev(center)[1],
                                        rev(df$limits.LCL)[1],
                                        rev(df$limits.UCL)[1]),
                                  x = rep(max(df$group)+0.5, 3)),
                aes_string(x = "x", y = "y"),
                label = c(label.center, label.limits),
                hjust = 0, nudge_x = 0.2,
                size = 10 * 5/14, col = gray(0.3))
  }
  
  # draw center line
  plot <- plot +
    geom_hline(yintercept = center, col = qcc.options("zones")$col[1])
  
  if(chart.all & (!is.null(newstats)))
  {
    len.obj.stats <- length(stats)
    len.new.stats <- length(newstats)
    plot <- plot +
      geom_vline(xintercept = min(xlim) + len.obj.stats + 0.5, lty = 3) +
      annotate("text", x = min(xlim) + len.obj.stats/2, 
               y = max(extendrange(ylim)),
               label = "Calibration data", 
               hjust = 0.5, vjust = -0.5, size = 10 * 5/14) +
      annotate("text", x = min(xlim) + len.obj.stats + len.new.stats/2,
               y = max(extendrange(ylim)),
               label = "New data", 
               hjust = 0.5, vjust = -0.5, size = 10 * 5/14)
  }
  
  if(add.stats) 
  { 
    # write info at bottom
    tab_base <- ggplot() + 
      ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + 
      theme_void() +
      theme(plot.background = element_rect(fill = qcc.options("bg.margin"),
                                           color = qcc.options("bg.margin")),
            plot.margin = margin(0.5, 0, 0.5, 0, unit = "lines"))

    text1 <- paste(paste0("Number of groups = ", length(statistics)),
                   paste0("Center = ", if(length(center) == 1) 
                     signif(center[1], digits) else "variable"),
                   paste0("StdDev = ", if(length(std.dev) == 1) 
                     signif(std.dev[1], digits) else "variable"), sep = "\n")
    
    text2 <- paste(paste0("Smoothing parameter = ", 
                          signif(object$lambda, digits = digits)),
                   paste0("Control limits at ", object$nsigmas, "xStdErr"),
                   paste0("No. beyond limits = ", 
                          sum(violations, na.rm = TRUE)), sep = "\n")
    tab1 <- tab_base + 
      geom_text(aes(x = -Inf, y = Inf), label = text1, 
                hjust = 0, vjust = 1, size = 10 * 5/14)
      # TODO: remove
      # theme(plot.margin = margin(0.5, 0, 0.5, 5, unit = "lines"))
    tab2 <- tab_base + 
      geom_text(aes(x = -Inf, y = Inf), label = text2, 
                hjust = 0, vjust = 1, size = 10 * 5/14)
      # TODO: remove
      # theme(plot.margin = margin(0.5, 1, 0.5, 3, unit = "lines"))
    # TODO: remove
    # plot <- gridExtra::arrangeGrob(plot, tab1, tab2,
    #                                layout_matrix = matrix(c(1,2,1,3), 
    #                                                       nrow = 2, ncol = 2),
    #                                heights = c(0.85, 0.15), 
    #                                widths = c(0.5, 0.5))
    plot <- patchwork::wrap_plots(plot, tab1, tab2,
                                  design = c("AA\nBC"),
                                  heights = c(0.85, 0.15), 
                                  widths = c(0.6, 0.4))
  }
  
  # class(plot) <- c("qccplot", class(plot))
  return(plot)
}

