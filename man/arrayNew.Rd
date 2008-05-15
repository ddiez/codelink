\name{arrayNew}
\alias{arrayNew}
\title{Create a new x11 device}
\description{
  Create a new x11 device with dimensions suited to be used
    with imageCodelink().
}
\usage{
  arrayNew(f=2, chip="rwgcod")
}
\arguments{
  \item{chip}{character; Codelink chip to be used in imageCodelink.}
  \item{f}{numerical; scaling factor.}
}
\author{Diego Diez}
\examples{
\dontrun{
    data(codelink.example)
    arrayNew()
    imageCodelink(codelink.example)
}
}
\keyword{documentation}
\keyword{utilities}
