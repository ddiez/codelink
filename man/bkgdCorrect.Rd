\name{bkgdCorrect}
\alias{bkgdCorrect}
\title{Background correction of intensity values.}
\description{
	Takes a Codelink object with Spot mean and Bkgd median values
	and corrects applying one of the methods available.
}
\usage{
  bkgdCorrect(object, method="half", preserve=FALSE,
              verbose=TRUE, offset=0)
}
\arguments{
  \item{object}{character; an object of class "Codelink".}
  \item{method}{character; the correction method to use, one of "none", "subtract", "half" and "normexp".}
  \item{preserve}{logical; if Smean and Bmedian slots should be preserved.}
  \item{verbose}{logical; if TRUE print some information with method normexp.}
  \item{offset}{numeric; value to add to intensities.}
}
\value{An object of class Codelink with corrected intensity values, that is
	Ri slot.}
\details{
  Available methods are:
  . None: Left intensities untouched.
  . Subtract: Simple subtraction of Bkgd median from Spot mean.
  . Half: The same as above but aboid negative values setting all intensity
          values below zero to 0.5.
  . Normexp: Apply normexp background adjustment from package limma.
}
\note{
}
\author{Diego Diez}
\seealso{
}
\examples{
#    data(codelink.example)
#    codelink.example <- bkgdCorrect(codelink.example, method="half")
}
\keyword{documentation}
\keyword{utilities}
