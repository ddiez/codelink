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
\details{
}
\note{
}
\author{Diego Díez Ruiz}
\seealso{\code{\link{log2}}}
\examples{
    #### NOT RUN
    ## data(codelink.example)
    ## codelink.example <- logCodelink(codelink.example)
}
\keyword{documentation}
\keyword{utilities}
