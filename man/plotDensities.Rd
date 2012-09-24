\name{plotDensities}
\alias{plotDensities}
\alias{plotDensities,Codelink-method}
\title{Plot Densities}
\description{
	Takes a Codelink object and plot the distributions of intensities.
}
\usage{
plotDensities(object, what = NULL, title = NULL, col = NULL, 
	legend.title = NULL, legend.cex=1, ...)
}
\arguments{
  \item{object}{an object of class "Codelink"}
  \item{what}{what data to plot, one of "bg", "smean", "snr", "ri" or "ni"}
  \item{title}{title of the plot}
  \item{col}{vector of colors}
  \item{legend.title}{if provided, a title for the legend}
  \item{legend.cex}{font factor use in legend}
  \item{...}{further arguments to be passed to plot() and line() functions (e.g. lwd, etc.)}
}
\author{Diego Diez}
\examples{
\dontrun{
    data(codelink.example)
    plotDensities(codelink.example)
}
}
\keyword{documentation}
\keyword{utilities}
