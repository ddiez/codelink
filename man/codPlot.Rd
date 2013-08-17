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

Plot types:

\describe{
    \item{\code{scatter}}{
    	Scatter plot of the selected array intensities vs. the median array, or a second array if specified.
    }
    \item{\code{ma}}{
    	MA plot of the selected array intensities vs. the median array, or a second array if specified. In the x-axis the value of A (mean abundance; A = Array1 + Array2) is plotted. In the y-axis the value of M (difference; M = Array1 - Array2) is plotted.
    }
    \item{\code{density}}{
    	The distribution of intensities for each array in the dataset is plotted as a kernel density.
    }
    \item{\code{image}}{
    	If the logical position of the probes in the array are available, this plot produces a pseudoimage of the array. The values of the "intensity", "bg" or "snr" can be selected with the 'signal' parameter (default signal="bg"). Color range can be controled with 'high' and 'low' parameters.
    }
	}
}
\author{
Diego Diez
}