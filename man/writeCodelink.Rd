\name{writeCodelink}
\alias{writeCodelink}
\title{Write a Codelink object to file.}
\description{
}
\usage{
  writeCodelink(object, file, dec = ".", sep = "\t", flag = FALSE)
}
\arguments{
  \item{object}{an object of class "Codelink".}
  \item{file}{filename to write object to.}
  \item{dec}{decimal character to use.}
  \item{sep}{delimiter character to use.}
  \item{flag}{should the Codelink flags be written.}
}
\value{
}
\details{
By default, intensities and SNR are wrote to the file. If set, the flag are 
also output. The header have  "INTENSITY\_", "SNR\_" and "FLAG\_" respectibely 
appendend to the sample name for those values. The default delimiter is the tab 
character, but that can be set with the sep argument. The default decimal 
character is the point.
}
\note{
}
\author{Diego Diez}
\seealso{
}
\examples{
\dontrun{
    data(codelink.example)
    writeCodelink(codelink.example, file = "foo.txt")
}
}
\keyword{documentation}
\keyword{utilities}
