\name{readHeader.Codelink}
\alias{readHeader.Codelink}
\title{Read Header from Codelink Bioarrays Files}
\description{
    Read the header of Codelink files and obtein usefull information.
}
\usage{
readHeader.Codelink(file, dec=FALSE)
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
\author{Diego Díez Ruiz}
\seealso{}
\examples{
    ## NOT RUN ##
    #files <- list.files(pattern="TXT")
    #head <- readHeader.Codelink(files[1])
    #############
}
\keyword{documentation}
\keyword{utilities}
