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
    \item{\code{Sample\_name}:}{Object of class \code{"character"} containing the probeset ID}
    \item{\code{pm}:}{Object of class \code{"matrix"} containing the PM
      intensities. Columns represent samples and rows the different probes.}
    \item{\code{mm}:}{Object of class \code{"matrix"} containing the MM intensities}
  }
}

\section{Methods}{
  \describe{
    \item{colnames}{\code{signature(x = "ProbeSet")}: the column names
      of the \code{pm} matrices which  are the sample names}
    \item{express.summary.stat}{\code{signature(x = "ProbeSet",
	pmcorrect = "character", summary = "character")}: applies a summary statistic to the probe set.}
    \item{sampleNames}{\code{signature(object = "ProbeSet")}: the column names
      of the \code{pm} matrices which  are the sample names}
  }
}
\note{More details are cotained in the vignette}

\seealso{}

\examples{
data(codelink.example)
}
\keyword{classes}
