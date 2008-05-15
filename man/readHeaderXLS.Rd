\name{readHeaderXLS}
\alias{readHeaderXLS}
\title{Read Header from XLS exported Codelink Bioarrays Files}
\description{
    Read the header of Codelink files and obtain useful information.
}
\usage{
readHeaderXLS(file, dec=FALSE)
}
\arguments{
  \item{file}{File to read.}
  \item{dec}{logical; If TRUE determine decimal point.}
}
\value{A list with header and other usefull information.}
\details{
	This function is not meant to be used by normal users.
}
\author{Diego Diez}
\examples{
\dontrun{
    files <- list.files(pattern = "TXT")
    head <- readHeaderHeader(files[1])
}
}
\keyword{documentation}
\keyword{utilities}
