\name{log2.Codelink}
\alias{log2.Codelink}
\title{Application of log2 to Codelink object}
\description{
	Takes a Codelink object and apply log2 to intensity values.
}
\usage{
log2.Codelink(object)
}
\arguments{
  \item{object}{an object of class "Codelink" or a list of genes.}
}
\value{A Codelink object with log2 intensities.}
\details{
}
\note{
}
\author{Diego Díez Ruiz}
\seealso{\code{\link{log2}}}
\examples{
    data(codelink.example)
    codelink.example <- log2.Codelink(codelink.example)
}
\keyword{documentation}
\keyword{utilities}
