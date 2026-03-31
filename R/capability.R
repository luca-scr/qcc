#-------------------------------------------------------------------#
#                                                                   #
#    Process Capability Analysis                                    #
#                                                                   #
#-------------------------------------------------------------------#

# TODO: implement Cpkm
processCapability <- function(object, spec.limits, target, 
                              std.dev, nsigmas, 
                              confidence.level = 0.95, ...)
{
# Computes process capability indices for a qcc object of type "xbar" 
# and plot the histogram

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")
  if (!(object$type=="xbar" | object$type=="xbar.one"))
     stop("Process Capability Analysis only available for charts type \"xbar\" and \"xbar.one\" charts")

  x <- as.vector(object$data)
  x <- x[!is.na(x)]
  center <- object$center
  if(missing(std.dev))
    std.dev <- object$std.dev
  overall.std.dev <- stats::sd(x)
  n <- length(x)

  if (missing(spec.limits))
     stop("specification limits must be provided")
  spec.limits <- as.vector(spec.limits)[1:2]
  LSL <- spec.limits[1]
  if(!(is.numeric(LSL) & is.finite(LSL))) LSL <- NA
  USL <- spec.limits[2]
  if(!(is.numeric(USL) & is.finite(USL))) USL <- NA
  if(is.na(LSL) & is.na(USL))
     stop("invalid specification limits")

  # TODO: refactor
  has.target <- (!missing(target))
  if(!has.target) {
    target <- mean(spec.limits, na.rm=TRUE) # FIX: removing NA means target == spec limits
    if(!is.na(LSL) & !is.na(USL)) has.target <- TRUE
    message("target value not provided; using midpoint of specification limits. Cpm and Ppm may be optimistic")
    # TODO: Explain this message in more detail in ?processCapability
  }
     
  if (is.na(LSL))
     { if (target > USL)
           warning("target value larger than one-sided specification limit...") }
  if (is.na(USL))
     { if (target < LSL)
           warning("target value smaller than one-sided specification limit...") }
  if (!is.na(LSL) & !is.na(USL))
     { if (target < LSL || target > USL)
       warning("target value is not within specification limits...") }
       
  if (missing(nsigmas))
     if (is.null(object$nsigmas))
        stop("nsigmas not available in the 'qcc' object. Please provide nsigmas.") 
     else  nsigmas <- object$nsigmas
  
  if (confidence.level < 0 | confidence.level > 1)
     stop("the argument confidence.level must be a value between 0 and 1") 

  # computes process capability indices
  Cp <- (USL - LSL) / (2*nsigmas*std.dev)
  Cp.u <- (USL-center)/(nsigmas*std.dev)
  Cp.l <- (center-LSL)/(nsigmas*std.dev)
  Cp.k <- min(Cp.u, Cp.l, na.rm = TRUE)
  # Cpm <- (USL - LSL) / (2*nsigmas*sqrt(sum((x-target)^2)/(n-1)))
  Cpm <- Cp / sqrt(1+((center-target)/std.dev)^2)
  Pp <- (USL - LSL) / (2*nsigmas*overall.std.dev)
  Pp.u <- (USL-center)/(nsigmas*overall.std.dev)
  Pp.l <- (center-LSL)/(nsigmas*overall.std.dev)
  Pp.k <- min(Pp.u, Pp.l, na.rm = TRUE)
  Ppm <- Pp / sqrt(1+((center-target)/overall.std.dev)^2)

  # compute confidence limits 
  alpha <- 1-confidence.level
  Cp.limits   <- .chisq_limits_cp_family(Cp, n - 1, alpha)
  Cp.u.limits <- .wald_limits_cpk_family(Cp.u, qnorm(confidence.level), n)
  Cp.l.limits <- .wald_limits_cpk_family(Cp.l, qnorm(confidence.level), n)
  Cp.k.limits <- .wald_limits_cpk_family(Cp.k, qnorm(1 - alpha / 2), n)
  df <- n * (1 + ((center - target) / std.dev)^2) / (1 + 2 * ((center - target) / std.dev)^2)
  Cpm.limits <- .chisq_limits_cp_family(Cpm, df, alpha)
  Pp.limits <- .chisq_limits_cp_family(Pp, n - 1, alpha)
  Pp.u.limits <- .wald_limits_cpk_family(Pp.u, qnorm(confidence.level), n)
  Pp.l.limits <- .wald_limits_cpk_family(Pp.l, qnorm(confidence.level), n)
  Pp.k.limits <- .wald_limits_cpk_family(Pp.k, qnorm(1 - alpha / 2), n)
  overall.df <- n * (1 + ((center - target) / overall.std.dev)^2) /
    (1 + 2 * ((center - target) / overall.std.dev)^2)
  Ppm.limits <- .chisq_limits_cp_family(Ppm, overall.df, alpha)

  limit.names <- c(paste(round(100*alpha/2, 1), "%", sep=""),
                   paste(round(100*(1-alpha/2), 1), "%", sep=""))
  names(Cp.limits) <- names(Cp.u.limits) <- names(Cp.l.limits) <- names(Cp.k.limits) <-
    names(Cpm.limits) <- names(Pp.limits) <- names(Pp.u.limits) <- names(Pp.l.limits) <-
    names(Pp.k.limits) <- names(Ppm.limits) <- limit.names

  if(is.na(LSL))  exp.LSL <- NA
  else { exp.LSL <- pnorm((LSL-center)/std.dev) * 100
         if(exp.LSL < 0.01) exp.LSL <- 0 }
  if(is.na(USL))  exp.USL <- NA
  else { exp.USL <- (1-pnorm((USL-center)/std.dev)) * 100
         if(exp.USL < 0.01) exp.USL <- 0 }
  obs.LSL <- sum(x<LSL)/n * 100
  obs.USL <- sum(x>USL)/n * 100
  
  tab <- cbind(c(Cp, Cp.l, Cp.u, Cp.k, Cpm, Pp, Pp.l, Pp.u, Pp.k, Ppm),
               rbind(Cp.limits, Cp.l.limits, Cp.u.limits, 
                     Cp.k.limits, Cpm.limits, Pp.limits, Pp.l.limits, Pp.u.limits,
                     Pp.k.limits, Ppm.limits))
  rownames(tab) <- c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm",
                     "Pp", "Pp_l", "Pp_u", "Pp_k", "Ppm")
  colnames(tab) <- c("Value", names(Cp.limits))

  out <- list(data = x, data.name = object$data.name,
              center = center, std.dev = std.dev, overall.std.dev = overall.std.dev,
              has.target = has.target, target = target, 
              spec.limits = { sl <- c(LSL, USL)
                              names(sl) <- c("LSL", "USL")
                              sl },
              indices = tab, 
              exp = { exp <- c(exp.LSL, exp.USL)
                      names(exp) <- c("Exp < LSL", "Exp > USL")
                      exp }, 
              obs = { obs <- c(obs.LSL, obs.USL)
                      names(obs) <- c("Obs < LSL", "Obs > USL")
                      obs },
              nobs = length(object$data) )
  class(out) <- "processCapability"
  return(out)
}

print.processCapability <- function(x, digits = getOption("digits"), ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  cat(cli::rule(left = crayon::bold("Process Capability Analysis"), 
                width = min(getOption("width"),50)), "\n\n")
  
  cat(paste(formatC("Number of obs = ", width=16),
            formatC(length(object$data), width=12, flag="-"),
            formatC("Target = ", width=10),
            ifelse(object$has.target, 
                   formatC(object$target, digits=digits, flag="-"), ""),
            "\n", sep=""))
  cat(paste(formatC("Center        = ", width=16),
            formatC(object$center, digits=digits, width=12, flag="-"),
            formatC("LSL    = ", width=10),
            ifelse(is.na(object$spec.limits[1]), "",
                   formatC(object$spec.limits[1], digits=digits, flag="-")),
            "\n", sep=""))
  cat(paste(formatC("StdDev        = ", width=16),
            formatC(object$std.dev, digits=digits, width=12, flag="-"),
            formatC("USL    = ", width=10),
            ifelse(is.na(object$spec.limits[2]), "",
                   formatC(object$spec.limits[2], digits=digits, flag="-")),
            "\n", sep=""))
  cat(paste(formatC("Overall SD    = ", width=16),
            formatC(object$overall.std.dev, digits=digits, width=12, flag="-"),
            "\n", sep=""))
            
  indices <- object$indices
  names(dimnames(indices)) <- c("Capability indices", "")
  print(indices, digits = 3, na.print="", print.gap=2)
  
  cat("\n")
  cat(paste("Exp<LSL", ifelse(is.na(object$exp[1]), "\t", 
                              paste(signif(object$exp[1], digits=2), "%\t", sep="")), 
            "Obs<LSL", ifelse(is.na(object$obs[1]), "", 
                              paste(signif(object$obs[1], digits=2), "%\n", sep=""))))
  cat(paste("Exp>USL", ifelse(is.na(object$exp[2]), "\t", 
                              paste(signif(object$exp[2], digits=2), "%\t", sep="")),
            "Obs>USL", ifelse(is.na(object$obs[2]), "", 
                              paste(signif(object$obs[2], digits=2), "%\n", sep=""))))
  
  invisible()
}

summary.processCapability <- function(object, ...) 
  print.processCapability(object, ...)

plot.processCapability <- function(x, 
                                   add.stats = qcc.options("add.stats"),
                                   breaks = nclass.hist, 
                                   fill = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.5), 
                                   color = "white",
                                   title, xlab,
                                   digits = getOption("digits"),
                                   ...)
{
# Computes the operating-characteristic curves for the S-chart with nsigmas
# limits. The values on the vertical axis give the probability of not detecting
# a change from sigma to c*sigma on the first sample following the change.

  object <- x  # Argh.  Really want to use 'object' anyway
   if ((missing(object)) | (!inherits(object, "processCapability")))
     stop("an object of class `processCapability' is required")

  nobs <- length(object$data)
  indices <- object$indices[, 1]
  Cp   <- indices["Cp"]
  Cp_l <- indices["Cp_l"]
  Cp_u <- indices["Cp_u"]
  Cp_k <- indices["Cp_k"]
  Cpm  <- indices["Cpm"]
  Pp   <- indices["Pp"]
  Pp_l <- indices["Pp_l"]
  Pp_u <- indices["Pp_u"]
  Pp_k <- indices["Pp_k"]
  Ppm  <- indices["Ppm"]
  if(is.function(breaks))
    breaks <- breaks(object$data)
  breaks <- as.integer(breaks)
  h <- hist(object$data, breaks = breaks, plot=FALSE)
  xlim <- range(c(h$breaks, object$spec.limits, object$target), na.rm = TRUE)
  xlim <- extendrange(r = xlim, f = 0.1)
  x  <- seq(min(xlim), max(xlim), length=250)
  dx <- dnorm(x, object$center, object$std.dev)
  ylim <- extendrange(c(h$density, dx))
  xlim <- range(c(h$breaks, x))
  
  if(missing(title))
    title <- "Process capability analysis"

  plot <- ggplot() +
    geom_histogram(data = data.frame(data = object$data),
                   aes(x = data, y = after_stat(density)),
                   stat = "bin", breaks = h$breaks,
                   fill = fill, color = color) +
    geom_line(data = data.frame(x, dx), 
              aes(x = x, y = dx)) +
    labs(title = title, subtitle = "", y = "",
         x = if(missing(xlab)) object$data.name else xlab) +
    coord_cartesian(xlim = xlim, ylim = ylim,
                    expand = FALSE, clip = "off") +
    theme_light() + 
    theme(plot.background = element_rect(fill = qcc.options("bg.margin"),
                                         color = qcc.options("bg.margin")),
          panel.background = element_rect(fill = qcc.options("bg.figure")),
          plot.title = element_text(face = "bold", size = 11),
          plot.margin = margin(5, 5, 5, 5),
          axis.text.y = element_text(angle = 90, 
                                     margin = margin(l = 5, r = 5),
                                     hjust = 0.5, vjust = 0.5))
  
  plot <- plot +
    geom_vline(xintercept = object$spec.limits, lty = 2) +
    annotate("text", 
             x = object$spec.limits[1], 
             y = ylim[2],
             label = "LSL", 
             hjust = 0.5, vjust = -0.5, size = 10 * 5/14) +
    annotate("text", 
             x = object$spec.limits[2], 
             y = ylim[2],
             label = "USL", 
             hjust = 0.5, vjust = -0.5, size = 10 * 5/14) 
  
  if(object$has.target)
  { 
    plot <- plot + 
      geom_vline(xintercept = object$target, lty = 3) +
      annotate("text", 
               x = object$target, 
               y = ylim[2],
               label = "Target", 
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

    text1 <- c(paste0("Number of obs = ", nobs),
               paste0("Center = ", signif(object$center, digits)),
               paste0("StdDev = ", signif(object$std.dev, digits)),
               paste0("Overall SD = ", signif(object$overall.std.dev, digits)))
    text1 <- paste(text1, collapse = "\n")
    tab1 <- tab_base + 
      geom_text(aes(x = -Inf, y = Inf), label = text1, 
                hjust = 0, vjust = 1, size = 10 * 5/14)
    # TODO: remove
    # theme(plot.margin = margin(0.5, 0, 0.5, 2, unit = "lines"))
    
    text2 <- paste(paste0("Target = ", if(object$has.target) signif(object$target, digits) else ""),
                   paste0("LSL = ", signif(object$spec.limits[1], digits)),
                   paste0("USL = ", signif(object$spec.limits[2], digits)), 
                   sep = "\n")
    tab2 <- tab_base + 
      geom_text(aes(x = -Inf, y = Inf), label = text2, 
                hjust = 0, vjust = 1, size = 10 * 5/14)
    # TODO: remove
    # theme(plot.margin = margin(0.5, 0, 0.5, 0.5, unit = "lines"))
    
    text3 <- paste(paste0("Cp     = ", ifelse(is.na(Cp), "", signif(Cp, 3))),
                   paste0("Cp_l  = ", ifelse(is.na(Cp_l), "", signif(Cp_l, 3))),
                   paste0("Cp_u = ", ifelse(is.na(Cp_u), "", signif(Cp_u, 3))),
                   paste0("Cp_k = ", ifelse(is.na(Cp_k), "", signif(Cp_k, 3))),
                   paste0("Cpm  = ", ifelse(is.na(Cpm), "", signif(Cpm, 3))),
                   sep="\n")
    tab3 <- tab_base + 
      geom_text(aes(x = -Inf, y = Inf), label = text3, 
                hjust = 0, vjust = 1, size = 10 * 5/14)
    # TODO: remove
    # theme(plot.margin = margin(0.5, 0, 0.5, 0.5, unit = "lines"))
    
    text4 <- paste(paste0("Pp     = ", ifelse(is.na(Pp), "", signif(Pp, 3))),
                    paste0("Pp_l  = ", ifelse(is.na(Pp_l), "", signif(Pp_l, 3))),
                    paste0("Pp_u = ", ifelse(is.na(Pp_u), "", signif(Pp_u, 3))),
                    paste0("Pp_k = ", ifelse(is.na(Pp_k), "", signif(Pp_k, 3))),
                    paste0("Ppm  = ", ifelse(is.na(Ppm), "", signif(Ppm, 3))),
                    sep="\n")
    tab4 <- tab_base + 
      geom_text(aes(x = -Inf, y = Inf), label = text4, 
                hjust = 0, vjust = 1, size = 10 * 5/14)

    text5 <- paste(paste0("Exp<LSL ", ifelse(is.na(object$exp[1]), "", paste0(signif(object$exp[1], 2), "%"))),
                   paste0("Exp>USL ", ifelse(is.na(object$exp[2]), "", paste0(signif(object$exp[2], 2), "%"))),
                   paste0("Obs<LSL ", ifelse(is.na(object$obs[1]), "", paste0(signif(object$obs[1], 2), "%"))),
                   paste0("Obs>USL ", ifelse(is.na(object$obs[2]), "", paste0(signif(object$obs[2], 2), "%"))),
                   sep="\n")
    tab5 <- tab_base + 
      geom_text(aes(x = -Inf, y = Inf), label = text5, 
                hjust = 0, vjust = 1, size = 10 * 5/14)
    # TODO: remove
    # theme(plot.margin = margin(0.5, 1, 0.2, 0.5, unit = "lines"))

    # TODO: remove
    # plot <- gridExtra::arrangeGrob(plot, tab1, tab2, tab3, tab4,
    #                                # gridExtra::grid.arrange(plot, tab1, tab2, tab3, tab4,
    #                                layout_matrix = matrix(c(1,2,1,3,1,4,1,5), 
    #                                                       nrow = 2, ncol = 4),
    #                                heights = c(0.78, 0.22), 
    #                                widths = c(0.35, 0.2, 0.2, 0.25))
    
    plot <- patchwork::wrap_plots(plotlist = list(plot, tab1, tab2, tab3, tab4, tab5),
                                  design = "AAAAA\nBCDEF",
                                  heights = c(0.73, 0.27), 
                                  widths = c(0.24, 0.16, 0.18, 0.18, 0.24))
  }

  # class(plot) <- c("qccplot", class(plot))
  return(plot)
}
  
# Compute approximate confidence limits for cpu, cpl, cpk.
# See @Bissell1990
# - idx: point estimate of the capability index
# - z:   normal quantile
# - n:   sample size
.wald_limits_cpk_family <- function(idx, z, n)
{
  if (is.na(idx)) return(c(NA_real_, NA_real_))
  idx * (1 + c(-1, 1) * z * sqrt(1 / (9 * n * idx^2) + 1 / (2 * (n - 1))))
}

# Compute confidence limits for cp and cpm.
# For Cp, see @Chou1990. For Cpm, see @Boyles1991
# - idx: point estimate of the capability index
# - df:  degrees of freedom
# - alpha: total tail probability
.chisq_limits_cp_family <- function(idx, df, alpha)
{
  if (is.na(idx)) return(c(NA_real_, NA_real_))
  idx * sqrt(qchisq(c(alpha / 2, 1 - alpha / 2), df) / df)
}
