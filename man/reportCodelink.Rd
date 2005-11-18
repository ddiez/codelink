\name{reportCodelink}
\alias{reportCodelink}
\title{Write a report of genes selected in HTML}
\description{
  Takes a list of genes as argument and writes an HTML page containing informa-
tion about these genes: Unigene, Genbank, Entrez Gene, etc.
}
\usage{
	reportCodelink(object, chip, filename = NULL, title = "Main title",
                 probe.type = FALSE, other = NULL, other.ord = NULL)
}
\arguments{
  \item{object}{an object of class "Codelink" or a list of genes.}
  \item{chip}{the chip description package.}
  \item{filename}{file name used in the report.}
  \item{title}{title used in the report.}
  \item{probe.type}{logical; if TRUE Probe type information is  writen.}
  \item{other}{list of vectors containing aditional values to add to the report.}
  \item{other.ord}{slot name in other to order genes by.}
}
\value{Nothing, only the HTML file generated.}
\details{
}
\note{
}
\author{Diego Díez Ruiz}
\seealso{\code{\link{htmlpage}}}
\examples{
}
\keyword{documentation}
\keyword{utilities}
