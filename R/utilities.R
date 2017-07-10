
blues.colors <- function (n) 
{
  palette <- grDevices::colorRampPalette(c("#03396c", "#005b96", "#6497b1", "#b3cde0"), 
                                         space = "Lab")
  palette(n)
}


#----------------------------------------------------------------------------#
# print a short version of a matrix by allowing to select the number of 
# head/tail rows and columns to display

.printShortMatrix <- function(x, head = 2, tail = 1, chead = 5, ctail = 1, ...)
{ 
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  if(is.na(head <- as.numeric(head))) head <- 2
  if(is.na(tail <- as.numeric(tail))) tail <- 1
  if(is.na(chead <- as.numeric(chead))) chead <- 5
  if(is.na(ctail <- as.numeric(ctail))) ctail <- 1
  
  if(nr > (head + tail + 1))
    { rnames <- rownames(x)
      if(is.null(rnames)) 
        rnames <- paste("[", 1:nr, ",]", sep ="")
      x <- rbind(x[1:head,,drop=FALSE], 
                 rep(NA, nc), 
                 x[(nr-tail+1):nr,,drop=FALSE])
      rownames(x) <- c(rnames[1:head], "...", rnames[(nr-tail+1):nr])
  }
  if(nc > (chead + ctail + 1))
    { cnames <- colnames(x)
      if(is.null(cnames)) 
        cnames <- paste("[,", 1:nc, "]", sep ="")
      x <- cbind(x[,1:chead,drop=FALSE], 
                 rep(NA, nrow(x)), 
                 x[,(nc-ctail+1):nc,drop=FALSE])
      colnames(x) <- c(cnames[1:chead], "...", cnames[(nc-ctail+1):nc])
  }
          
  print(x, na.print = "", ...)
}

# old version
# .printShortMatrix <- function(x, head = 2, tail = 1, ...)
# { 
#   x <- as.matrix(x)
#   nr <- nrow(x)
#   nc <- ncol(x)
#   if(nr > 4)
#     { rnames <- rownames(x)
#       if(is.null(rnames)) 
#         rnames <- paste("[", 1:nr, ",]", sep ="")
#       x <- rbind(x[1:head,], rep(NA, nc), x[(nr-tail+1):nr,])
#       rownames(x) <- c(rnames[1:head], "...", rnames[(nr-tail+1):nr])
#       print(x, na.print = "", ...)
#     }
#   else
#     { print(x, ...)}  
# }

#-------------------------------------------------------------------#
#                                                                   #
#
# Options retrieval and setting
#

qcc.options <- function (...)
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

".qcc.options" <- list(exp.R.unscaled = c(NA, 1.128, 1.693, 2.059, 2.326, 2.534, 2.704, 2.847, 2.970, 3.078, 3.173, 3.258, 3.336, 3.407, 3.472, 3.532, 3.588, 3.640, 3.689, 3.735, 3.778, 3.819, 3.858, 3.895, 3.931),
                       se.R.unscaled = c(NA, 0.8525033, 0.8883697, 0.8798108, 0.8640855, 0.8480442, 0.8332108, 0.8198378, 0.8078413, 0.7970584, 0.7873230, 0.7784873, 0.7704257, 0.7630330, 0.7562217, 0.7499188, 0.7440627, 0.7386021, 0.7334929, 0.7286980, 0.7241851, 0.7199267, 0.7158987, 0.7120802, 0.7084528, 0.7050004, 0.7017086, 0.6985648, 0.6955576, 0.6926770, 0.6899137, 0.6872596, 0.6847074, 0.6822502, 0.6798821, 0.6775973, 0.6753910, 0.6732584, 0.6711952, 0.6691976, 0.6672619, 0.6653848, 0.6635632, 0.6617943, 0.6600754, 0.6584041, 0.6567780, 0.6551950, 0.6536532, 0.6521506),
                      beyond.limits = list(pch=19, col="red"),
                      violating.runs = list(pch=19, col="orange"),
                      run.length = 7,
                      # bg.margin = "lightgrey",
                      bg.margin = "#E5E5E5",
                      bg.figure = "white",
                      cex = 1,
                      font.stats = 1,
                      cex.stats = 0.9)

.onAttach <- function(library, pkg)
{
  ## we can't do this in .onLoad
  unlockBinding(".qcc.options", asNamespace("qcc"))
  description <- readLines(system.file("DESCRIPTION", package = "qcc"))
  version <- grep("Version:", description, ignore.case = TRUE, value = TRUE)
  version <- gsub(pattern = "Version:", replacement = "", version, ignore.case = TRUE)
  version <- gsub(pattern = " ", replacement = "", version)
  packageStartupMessage("Package 'qcc', version ", version)
  packageStartupMessage("Type 'citation(\"qcc\")' for citing this R package in publications.")
  invisible()
}

