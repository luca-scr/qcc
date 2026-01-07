#-------------------------------------------------------------------#
#                                                                   #
#                      CUSUM CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

cusum <- function(data, 
                  sizes, center, std.dev, 
                  decision.interval = 5, se.shift = 1,
                  head.start = 0, 
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

  if (decision.interval <= 0)
      stop("decision.interval must be positive")

  if (head.start < 0 || head.start >= decision.interval)
      stop("head.start must be non-negative and less than decision.interval")
  
  # used for computing statistics and std.dev
  type <- if(any(sizes==1)) "xbar.one" else "xbar"
  if(ncol(data) == 1 & any(sizes > 1) & missing(std.dev))
     stop("sizes larger than 1 but data appears to be single samples. In this case you must provide also the std.dev")
  
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
  names(statistics) <-  rownames(data) <-  labels
  names(dimnames(data)) <- list("Group", "Samples")

  object <- list(call = call, type = "cusum", 
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
  violations <- list(lower = ifelse(cusum.neg < ldb, 1, NA),
                     upper = ifelse(cusum.pos > udb, 1, NA))
  
  object$type <- "cusum"
  object$pos <- cusum.pos 
  object$neg <- cusum.neg
  object$head.start <- head.start
  object$decision.interval <- decision.interval
  object$se.shift <- se.shift
  object$violations <- violations

  class(object) <- "cusum.qcc"
  return(object)
}

print.cusum.qcc <- function(x, digits =  getOption("digits"), ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  # cat("\nCall:\n",deparse(object$call),"\n\n",sep="")
  cat(cli::rule(left = crayon::bold("Cusum Chart"), 
                width = min(getOption("width"),50)), "\n\n")
  
  data.name <- object$data.name
  # type <- object$type
  statistics <- object$statistics
  # cat(paste(type, "chart for", data.name, "\n"))
  # cat("\nSummary of group statistics:\n")
  # print(summary(statistics), ...)
  
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
  
  cat("\n")
  if(object$head.start > 0)
  {
    cat("Head start (StdErr)        =",
        signif(object$head.start, digits = digits), "\n")
  }
  cat("Decision interval (StdErr) =", 
      signif(object$decision.interval, digits = digits), "\n")
  cat("Shift detection   (StdErr) =", 
      signif(object$se.shift, digits = digits), "\n")

  invisible()
}

summary.cusum.qcc <- function(object, ...) print.cusum.qcc(object, ...)

plot.cusum.qcc <- function(x, xtime = NULL,
                           add.stats = qcc.options("add.stats"), 
                           chart.all = qcc.options("chart.all"), 
                           fill = qcc.options("fill"),
                           label.bounds = c("LDB", "UDB"), 
                           title, xlab, ylab, xlim, ylim,
                           digits = getOption("digits"), 
                           ...) 
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
                   cusum_pos = cusum.pos,
                   cusum_neg = cusum.neg,
                   ldb = ldb, udb = udb,
                   violations_lower = factor(ifelse(is.na(violations$lower), 
                                                    0, violations$lower),
                                             levels = 0:1),
                   violations_upper = factor(ifelse(is.na(violations$upper), 
                                                    0, violations$upper),
                                             levels = 0:1))
  if(!chart.all & (!is.null(newstats)))
    df <- df[df$group > length(object$statistics),]

  if(missing(ylim))
    ylim <- range(df[,c("cusum_pos", "cusum_neg", "ldb", "udb")], na.rm = TRUE)
  if(missing(xlim))
    xlim <- range(df$group, na.rm = TRUE)
  
  plot <- 
    ggplot(data = df, aes(x = .data[["group"]])) +
    geom_line(aes(y = .data[["cusum_pos"]])) +
    geom_point(aes(y = .data[["cusum_pos"]], 
                   colour = .data[["violations_upper"]], 
                   shape = .data[["violations_upper"]]), 
               size = 2) +
    geom_line(aes(y = .data[["cusum_neg"]])) +
    geom_point(aes(y = .data[["cusum_neg"]], 
                   colour = .data[["violations_lower"]], 
                   shape = .data[["violations_lower"]]), 
               size = 2) +
    scale_colour_manual(values = c("black", qcc.options("rules")$col),
                        breaks = levels(df$violations)) +
    scale_shape_manual(values = c(20, qcc.options("rules")$pch),
                       breaks = levels(df$violations)) +
    labs(title = title, subtitle = "",
         x = if(missing(xlab)) "Group" else xlab,
         y = if(missing(ylab)) "Cumulative Sum" else ylab) +
    coord_cartesian(xlim = xlim+c(-0.5,0.5), 
                    ylim = extendrange(ylim),
                    expand = FALSE, clip = "off") +
    theme_light() + 
    theme(plot.background = element_rect(fill = qcc.options("bg.margin"),
                                         color = qcc.options("bg.margin")),
          panel.background = element_rect(fill = qcc.options("bg.figure")),
          plot.title = element_text(face = "bold", size = 11),
          axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
          legend.position = "none",
          plot.margin = margin(5, 30, 5, 5),
          axis.text.y = element_text(angle = 90, 
                                     margin = margin(l = 5, r = 5),
                                     hjust = 0.5, vjust = 0.5))
  
  plot <- plot + 
  {
    if(is.numeric(groups))
      scale_x_continuous(breaks = pretty(df$group, n = 7))
    else
      scale_x_date(breaks = pretty(df$group, n = 7))
  }
    
  lab <- "Above target"
  if (add.stats && object$head.start > 0)
      lab <- paste(lab, " (start = ", object$head.start, ")", sep = "")
  plot <- plot + 
    annotate("text", x = min(xlim)-0.1*diff(range(xlim)),
             y = max(extendrange(ylim))/2,
             label = lab, angle = 90,
             col = gray(0.3), size = 10 * 5/14,
             hjust = 0.5, vjust = 0.5)
  lab <- "Below target"
  if (add.stats && object$head.start > 0)
    lab <- paste(lab, " (start = ", - object$head.start, ")", sep = "")
  plot <- plot + 
    annotate("text", x = min(xlim)-0.1*diff(range(xlim)),
             y = min(extendrange(ylim))/2,
             label = lab, angle = 90,
             col = gray(0.3), size = 10 * 5/14,
             hjust = 0.5, vjust = 0.5)
    
  # draw decision boundaries
  if(all(is.finite(ldb)) & is.finite(udb))
  { 
    if(fill)
    { 
      xp <- c(min(df$group)-0.5, max(df$group)+0.5)
      xp <- c(xp, rev(xp))
      yp <- c(ldb, ldb, udb, udb)
      
      plot <- plot + 
        geom_polygon(data = data.frame(xp, yp),
                     aes(x = .data[["xp"]], 
                         y = .data[["yp"]]), 
                     fill = adjustcolor(qcc.options("zones")$fill, alpha.f=0.2),
                     col = NA)
    } else
    {
      plot <- plot + 
        geom_hline(yintercept = ldb,
                  lty = qcc.options("zones")$lty[1],
                  col = qcc.options("zones")$col[1])
      plot <- plot + 
        geom_hline(yintercept = udb, 
                  lty = qcc.options("zones")$lty[1],
                  col = qcc.options("zones")$col[1])
    }
  }
  
  # draw center line
  plot <- plot +
    geom_hline(yintercept = 0, col = qcc.options("zones")$col[1])

  if(chart.all & (!is.null(newstats)))
  {
    len.obj.stats <- length(stats)
    len.new.stats <- length(newstats)
    plot <- plot +
      geom_vline(xintercept = mean(groups[len.obj.stats+c(0,1)]), lty = 3) +
      # geom_vline(xintercept = min(xlim) + len.obj.stats + 0.5, lty = 3) +
      annotate("text", 
               # x = min(xlim) + len.obj.stats/2, 
               x = mean(c(groups[1], mean(groups[len.obj.stats+c(0,1)]))),
               # y = max(extendrange(ylim)),
               y = max(ylim),
               label = "Calibration data", 
               hjust = 0.5, vjust = -0.5, size = 10 * 5/14) +
      annotate("text", 
               # x = min(xlim) + len.obj.stats + len.new.stats/2,
               x = mean(c(mean(groups[len.obj.stats+c(0,1)]), groups[len.obj.stats+len.new.stats])),
               # y = max(extendrange(ylim)),
               y = max(ylim),
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
    
    text2 <- paste(paste0("Decision interval (StdErr) = ", 
                          signif(object$decision.interval, digits = digits)),
                   paste0("Shift detection (StdErr) = ", 
                          signif(object$se.shift, digits = digits)),
                   paste0("No. beyond boundaries = ", 
                          sum(sapply(violations, sum, na.rm = TRUE))), sep = "\n")
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

