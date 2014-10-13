\name{CodelinkSetUnique-class}
\alias{class::CodelinkSetUnique}
\alias{CodelinkSetUnique}
\alias{CodelinkSetUnique-class}
\docType{class}

\title{Class CodelinkSetUnique}
\description{This is the storage class for Codelink data-- with unique probes}
\section{Objects from the Class}{
  Object are created after applying averageProbes() on an object of the class CodelinkSet
}

\section{Description}{
	CodelinkSetUnique objects are derived from ExpressionSet and therefore inherits all the methods.}

\note{More details are in the package vignette}

\author{Diego Diez}

\examples{
\dontrun{
    data(codelink.exprs)
    foo <- averageProbes(codelink.exprs)
}
}
\keyword{classes}
