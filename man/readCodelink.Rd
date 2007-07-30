\name{readCodelink}
\alias{readCodelink}
\title{Read Codelink Bioarrays Data}
\description{
  Read data exported as text by Codelink Software. It reads values (normalized
by Codelink Software or not) flags and information about probes.
}
\usage{
    readCodelink(files=list.files(pattern="TXT"), sample.name=NULL,
	flag=list(M=NA, I=NA, C=NA), dec=NULL, type="Spot",	preserve=FALSE,
	verbose=2, file.type="Codelink", check=TRUE, fix=FALSE)
}
\arguments{
  \item{files}{list of files to read.}
  \item{sample.name}{vector of same length as files with sample names.}
  \item{flag}{list with values to assign based on Flag quality values.}
  \item{dec}{character indicating the decimal character used in the files.}
  \item{type}{character indicating which base value to read from files.}
  \item{preserve}{logical, if TRUE Bkgd\_stdev slot is not removed (if present).}
  \item{verbose}{numerical, set the level of information. Level 2 set as old behaviour. Level > 2 output some debug info.}
  \item{file.type}{exported file type, currently Codelink or XLS file formats supported.}
  \item{check}{logical, check for probe order consistency.}
  \item{fix}{logical, try to fix probe order consistency.}
}
\value{An object of class "Codelink".}
\details{
}
\note{
}
\author{Diego Diez}
\seealso{\code{\link{read.table}}}
\examples{
\dontrun{
    # reading default extension (TXT).
    data <- readCodelink()

    # specify a different one.
    files <- list.files(pattern = "txt")
    data <- readCodelink(files = files)

    # example.
    data(codelink.example)
}
}
\keyword{documentation}
\keyword{utilities}
