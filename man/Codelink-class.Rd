\name{Codelink-class}
\docType{class}
\alias{class:Codelink}
\alias{Codelink}
\alias{Codelink-class}
\alias{[,Codelink-method}
\alias{[,Codelink,ANY,ANY-method}
\alias{show,Codelink-method}

\title{Class Codelink}
\description{This is the storage class for Codelink data}
\section{Objects from the Class}{
  Object are created after reading text codelink files with readCodelink()}

\section{Description}{
	Codelink objects contain a single \code{"list"}, which contains the following elements:
  \describe{
    \item{sample}{Object of class \code{"character"} containing the sample names}
    \item{file}{Object of class \code{"character"} containing the file names}
	\item{name}{Object of class \code{"character"} containing the probe ID}
    \item{method}{Object of class \code{"list"} containing log information}
    \item{Smean}{Object of class \code{"matrix"} containing spot mean intensities}
    \item{Bmedian}{Object of class \code{"matrix"} containing background median intensities}
    \item{Ri}{Object of class \code{"matrix"} containing raw intensities}
    \item{Ni}{Object of class \code{"matrix"} containing normalized intensities}
    \item{snr}{Object of class \code{"matrix"} containing signal to noise ratio values}
    \item{flag}{Object of class \code{"character"} containing assigned flags}
  }
}

\note{More details are in the package vignette}

\author{Diego Diez}

\examples{
\dontrun{
    data(codelink.example)
}
}
\keyword{classes}
