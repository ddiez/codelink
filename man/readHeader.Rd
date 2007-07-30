\name{readHeader}
\alias{readHeader}
\title{Read Header from Codelink Bioarrays Files}
\description{
    Read the header of Codelink files and obtein usefull information.
}
\usage{
readHeader(file, dec=FALSE)
}
\arguments{
  \item{file}{File to read.}
  \item{dec}{logical; If TRUE determine decimal point.}
}
\value{A list with header and other usefull information.}
\details{
}
\note{
}
\author{Diego Diez}
\seealso{}
\examples{
\dontrun{
    files <- list.files(pattern = "TXT")
    head <- readHeader(files[1])
}
}
\keyword{documentation}
\keyword{utilities}
