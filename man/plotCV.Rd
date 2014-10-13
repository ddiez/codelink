\name{plotCV}
\alias{plotCV}
\title{Plot of CV}
\description{
	Takes a Codelink object and plot de distribution of CV after applying 
	mergeCodelink.
}
\usage{
	plotCV(object, subset=c(1:dim(object)[2]), cutoff=NULL, title=NULL, legend.cex=1)
}
\arguments{
	\item{object}{an object of class "Codelink".}
	\item{subset}{subset of arrays to plot}
	\item{cutoff}{cutoff of CV to be shown.}
	\item{title}{title of the plot.}
	\item{legend.cex}{factor to apply to the fonts in the legend to fit.}
}
\author{Diego Diez}
\keyword{documentation}
\keyword{utilities}
