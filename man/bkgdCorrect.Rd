\name{bkgdCorrect}
\alias{bkgdCorrect}
\alias{codCorrect}
\alias{codCorrect,CodelinkSet-method}

\title{Background correction of intensity values.}
\description{
Takes a CodelinkSet or Codelink object with Spot mean and Bkgd median values and performes background correction using one of the methods available.
}
\usage{
\method{bkgdCorrect}{Codelink}(object, method = "half", preserve = FALSE, verbose = FALSE, offset = 0)

\S4method{codCorrect}{CodelinkSet}(object, method = "half", normexp.method="saddle", offset=0)
}
\arguments{
  \item{object}{an object of class CodelinkSet or Codelink.}
  \item{method}{the correction method to use, one of "none", "subtract", "half" and "normexp".}
  \item{preserve}{logical; if Smean and Bmedian slots should be preserved.}
  \item{verbose}{logical; if TRUE print some information with method normexp.}
  \item{normexp.method}{method used with normexp correction.}
  \item{offset}{numeric; value to add to intensities.}
}
\value{An object of class Codelink with corrected intensity values, that is
	Ri slot.}
\details{
  Available methods are:
  . none: left intensities untouched.
  . subtract: simple subtraction of Bkgd median from Spot mean.
  . half: the same as above but aboid negative values setting all intensity
          values below zero to 0.5.
  . normexp: apply normexp background adjustment from package limma.
}
\author{Diego Diez}
\examples{
    data(codelink.exprset)
    codelink.exprset <- bkgdCorrect(codelink.exprset, method = "half")
}
\keyword{documentation}
\keyword{utilities}
