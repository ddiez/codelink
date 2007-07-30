\name{Codelink-class}
\docType{class}
\alias{Codelink-class}
\alias{[,Codelink-method}
\alias{show,Codelink-method}

\title{Class Codelink}
\description{This is the storage class for Codelink data}
\section{Objects from the Class}{
  Object are created after reading text codelink files with read.Codelink()}
\section{Slots}{
  \describe{
    \item{\code{sample}:}{Object of class \code{"character"} containing the sample names}
    \item{\code{file}:}{Object of class \code{"character"} containing the file names}
	\item{\code{name}:}{Object of class \code{"character"} containing the probe ID}
    \item{\code{method}:}{Object of class \code{"list"} containing log information}
    \item{\code{Smean}:}{Object of class \code{"matrix"} containing spot mean intensities}
    \item{\code{Bmedian}:}{Object of class \code{"matrix"} containing background median intensities}
    \item{\code{Ri}:}{Object of class \code{"matrix"} containing raw intensities}
    \item{\code{Ni}:}{Object of class \code{"matrix"} containing normalized intensities}
    \item{\code{snr}:}{Object of class \code{"matrix"} containing signal to noise ratio values}
    \item{\code{flag}:}{Object of class \code{"character"} containing assigned flags}
  }
}

\section{Methods}{
  \describe{
  }
}
\note{More details are cotained in the vignette}

\author{Diego Diez}
\seealso{}

\examples{
\dontrun{
    data(codelink.example)
}
}
\keyword{classes}
