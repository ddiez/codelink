\name{bkgdCorrect}
\alias{bkgdCorrect}
\title{Background correction of intensity values.}
\description{
	Takes a Codelink object with Spot mean and Bkgd median values
	and corrects applying one of the methods available.
}
\usage{
	bkgdCorrect(object,method="half")
}
\arguments{
  \item{object}{an object of class "Codelink".}
  \item{method}{the correction method to use.}
}
\value{An object of class Codelink with corrected intensity values, that is
	Raw intensity slot.}
\details{
  Available methods are:
  . Subtract: Simple subtraction of Bkgd median from Spot mean.
  . Half: The same as above but aboid negative values setting all intensity
          values below zero to 0.5.
}
\note{
}
\author{Diego Diez Ruiz}
\seealso{
}
\examples{
    data(codelink.example)
    codelink.example <- bkgdCorrect(codelink.example, method="half")
}
\keyword{documentation}
\keyword{utilities}
