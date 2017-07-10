\name{cause.and.effect}
\alias{cause.and.effect}
\title{Cause and effect diagram}
\description{Draw a basic cause and effect diagram.}
\usage{
cause.and.effect(cause, effect, title = "Cause-and-Effect diagram", 
                 cex = c(1, 0.9, 1), font = c(1, 3, 2))
}
\arguments{
\item{cause}{a list of causes and branches providing descriptive labels (see the example below).}
\item{effect}{a string label or the effect.}
\item{title}{a string specifying the main title to appear on the plot.}
\item{cex}{a vector of values for the graphical character expansion. The values refer, in order, to branches, causes and effect.}
\item{font}{a vector of values for the font to use. The values refer, in order, to branches, causes and effect.}
}
%\details{}
%\value{}
\references{
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
\examples{
cause.and.effect(cause=list(Measurements=c("Micrometers", "Microscopes", "Inspectors"),
                            Materials=c("Alloys", "Lubricants", "Suppliers"),
                            Personnel=c("Shifts", "Supervisors", "Training", "Operators"),
                            Environment=c("Condensation", "Moisture"),
                            Methods=c("Brake", "Engager", "Angle"),
                            Machines=c("Speed", "Lathes", "Bits", "Sockets")),
                 effect="Surface Flaws")
}
\keyword{hplot}
