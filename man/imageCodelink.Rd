\name{imageCodelink}
\alias{imageCodelink}
\title{Image plot of Codelink arrays}
\description{
  Plot and image of a Codelink array if the layout information is found.
}
\usage{

imageCodelink(object, array = 1, what = "bg",
    low="black", high="white", mar=c(2,1,1,1),
    gr=1, gc=1, log.it=FALSE, ...)
	
}
\arguments{
  \item{object}{an object of class "Codelink".}
  \item{array}{array to be used.}
  \item{what}{with data plot: bg, smean, ri, ni.}
  \item{low}{color used for low intensities.}
  \item{high}{color used for high intensities-}
  \item{mar}{character vector specifiying margings.}
  \item{gc}{numerical; number of grid columns.}
  \item{gr}{numerical; number of grid rows.}
  \item{log.it}{logical; if TRUE data is log2 transformed (if not yet).}
  \item{...}{additional arguments passed to image.}
}
\value{
}
\details{
}
\note{
}
\author{Diego Diez}
\seealso{
}
\examples{
#    data(codelink.example)
#    imageCodelink(codelink.example)
}
\keyword{documentation}
\keyword{utilities}
