\name{logCodelink}
\alias{logCodelink}
\title{Application of logCodelink to Codelink object}
\description{
	Takes a Codelink object and apply logCodelink to intensity values.
}
\usage{
logCodelink(object)
}
\arguments{
  \item{object}{an object of class "Codelink" or a list of genes.}
}
\value{A Codelink object with logCodelink intensities.}
\author{Diego Diez}
\seealso{\code{\link{log2}}}
\examples{
\dontrun{
    data(codelink.example)
    codelink.example <- logCodelink(codelink.example)
}
}
\keyword{documentation}
\keyword{utilities}
