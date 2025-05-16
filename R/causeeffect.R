#-------------------------------------------------------------------#
#                                                                   #
#                  Cause-and-Effect Diagram                         #
#                                                                   #
#-------------------------------------------------------------------#

causeEffectDiagram <- function(cause, effect, 
                               title = "Cause-and-Effect diagram",
                               cex = c(0.9,1,1.2), 
                               font = c(3,1,2),
                               ...)
{
  
  # running mean of successive pairs of obs
  mean2 <- function(x)
  { 
    m <- rep(NA, length(x)-1)
    for (i in 1:(length(x)-1))
      m[i] <- mean(x[c(i,i+1)])
    return(m)
  }
  
  nc <- length(cause)
  ncup <- nc - round(nc/2)
  nclo <- nc - ncup
  ncc <- max(sapply(cause, length))

  plot <- ggplot() +
    xlim(0, 100) + ylim(0, 100) +
    labs(title = title) +
    theme_void() +
    theme(plot.background = element_rect(fill = qcc.options("bg.margin"),
                                         color = qcc.options("bg.margin")),
          panel.background = element_rect(fill = qcc.options("bg.figure")),
          plot.title = element_text(face = "bold", margin = margin(b = 10)),
          plot.margin = margin(10, 10, 10, 10))
  
  size <- cex*12/ggplot2::.pt
  inches_per_unit <- 6 / 100  # ~ 6 inches across 100 units
  we <- strwidth(effect, units = "inches") /inches_per_unit * 1.1
  wc <- max(unlist(sapply(cause, strwidth, units = "inches"))) /inches_per_unit
  hc <- max(strheight(effect, units="inches"),
            unlist(sapply(cause, strheight, units="inches"))) /inches_per_unit

  plot <- plot +  
    # add main spine
    geom_segment(aes(x = 0, y = 50, xend = 100-we-2, yend = 50),
                 arrow = arrow(length = unit(0.02, "npc")), 
                 linewidth = 1) +
    # add effect label
    annotate("text", x = 100-we, y = 50,  
             label = effect, hjust = 0, 
             fontface = font[3], size = size[3])

  # draw branches and cause labels
  a <- (100-we)/(max(ncup,nclo)+1)
  ac <- a*(0:(max(ncup,nclo)))
  for(i in 1:(length(ac)-1))
  { 
    plot <- plot +
      annotate("segment", x = mean2(ac)[i], y = 95, xend = ac[i+1], yend = 50) +
      annotate("text", x = mean2(ac)[i], y = 97, label = names(cause)[i], 
               vjust = 0, fontface = font[2], size = size[2])
    if(i <= nclo)
    { 
      plot <- plot + 
        annotate("segment", x = mean2(ac)[i], y = 5, xend = ac[i+1], yend = 50) +
        annotate("text", x = mean2(ac)[i], y = 3, label = names(cause)[ncup+i], 
                 vjust = 1, fontface = font[2], size = size[2])
    }
  }
  # draw labels for upper branches
  for (j in 1:ncup)
  { 
    b <- (50-95)/(ac[j+1]-mean2(ac)[j])
    a <- 95-b*mean2(ac)[j]
    y <- rev(50+cumsum((95-50)/(ncc+1))*(1:(ncc)))
    x <- (y-a)/b
    for (i in 1:length(y))
    { 
      label <- cause[[j]][i]
      if (!is.na(label))
      {
        plot <- plot + 
          annotate("text", x = x[i], y = y[i], label = label, 
                   hjust = -0.1, fontface = font[1], size = size[1])
      }
    }
  }
  # draw labels for lower branches
  for (j in 1:ncup)
  { 
    b <- (50-5)/(ac[j+1]-mean2(ac)[j])
    a <- 5-b*mean2(ac)[j]
    y <- cumsum((95-50)/(ncc+1))*(1:(ncc))
    x <- (y-a)/b
    if (j <= nclo)
    {
      for (i in 1:length(y))
      { 
        label <- cause[[ncup+j]][i]
        if (!is.na(label))
        {
          plot <- plot + 
            annotate("text", x = x[i], y = y[i], label = label, 
                     hjust = -0.1, fontface = font[1], size = size[1])
        }  
      }
    }
  }
  
  return(plot)
}

