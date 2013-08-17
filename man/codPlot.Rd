\name{codPlot}
\alias{codPlot}
\alias{codPlot,CodelinkSet}
\alias{codPlot,MArrayLM-method}

\title{
Diagnostic plots for CodelinkSet object.
}
\description{
This function performs several types of diagnostic plots using information from CodelinkSet objects.
}
\usage{
codPlot(x, array, what = "ma", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a CodelinkSet object.
}
  \item{array}{
The array to use for plotting (by default array=1).
}
  \item{what}{
the type of plot (by default ma).
}
  \item{\dots}{
additional arguments passed to the appropriate functions.
}
}
\details{
ma - MA plot.
density - density plot.
scatter - scatter correlation plot.
image - array pseudo images.
}
\author{
Diego Diez
}