\name{CodelinkSet-class}
\alias{class::CodelinkSet}
\alias{CodelinkSet}
\alias{CodelinkSet-class}
\alias{codPreprocess}
\alias{codPreprocess,CodelinkSet-method}
\alias{getBkg}
\alias{getInt}
\alias{getInfo}
\alias{getSNR}
\alias{getFlag}
\alias{getWeight}
\alias{getBkg,CodelinkSet-method}
\alias{getInt,CodelinkSet-method}
\alias{getInfo,CodelinkSet-method}
\alias{getSNR,CodelinkSet-method}
\alias{getFlag,CodelinkSet-method}
\alias{getWeight,CodelinkSet-method}
\alias{exprs,CodelinkSet-method}
\alias{Codelink2CodelinkSet}


\alias{probeNames}
\alias{probeNames,CodelinkSet-method}
\alias{probeNames,CodelinkSetUnique-method}
\alias{probeTypes}
\alias{probeTypes,CodelinkSet-method}

\docType{class}

\title{Class CodelinkSet}
\description{This is the storage class for Codelink data}
\section{Objects from the Class}{
  Object are created after reading text codelink files with readCodelinkSet()}

\section{Description}{
	CodelinkSet objects are derived from ExpressionSet and therefore inherits all their methods. Additional methods are defined to extract information:

	Old Codelink-class objects can be converted into CodelinkSet instances using the function Codelink2CodelinkSet.
}

\section{Extends}{
   Directly extends class \code{\link{ExpressionSet}}.
}

\section{Methods}{
Class-specific methods:

	\describe{
    \item{\code{getWeight}}{Get the matrix of weights.}
	}
}

\author{Diego Diez}

\examples{
\dontrun{
    data(codset)
}
}
\keyword{classes}
