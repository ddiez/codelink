\name{readCodelink}
\alias{readCodelink}
\alias{readCodelinkSet}
\title{Read Codelink Bioarrays Data}
\description{
  Read data exported as text by Codelink Software. It reads values (normalized
by Codelink Software or not) flags and information about probes.
}
\usage{
readCodelink(files = list.files(pattern = "TXT"), sample.name =
                 NULL, flag, flag.weights, type.weights, dec = NULL,
                 type = "Spot", preserve = FALSE, verbose = 2,
                 file.type = "Codelink", check = TRUE, fix = FALSE, old
                 = FALSE)
readCodelinkSet(filename, columns = list(Signal = "Spot_mean", Background = "Bkgd_median"), phenoData = NULL, \dots)
}
\arguments{
  \item{files}{list of files to read.}
  \item{sample.name}{vector of same length as files with sample names.}
  \item{flag}{list with values to assign based on Flag quality values.}
  \item{flag.weights}{weights assigned to each probe flag.}
  \item{type.weights}{weights assigned to each probe type.}
  \item{dec}{character indicating the decimal character used in the files.}
  \item{type}{character indicating which base value to read from files.}
  \item{preserve}{logical, if TRUE Bkgd\_stdev slot is not removed (if present).}
  \item{verbose}{numerical, set the level of information. Level 2 set as old behaviour. Level > 2 output some debug info.}
  \item{file.type}{exported file type, currently Codelink or XLS file formats supported.}
  \item{check}{logical, check for probe order consistency.}
  \item{fix}{logical, try to fix probe order consistency.}
  \item{old}{logical, whether NA is assigned to intensities based on flags.}
  
  \item{filename}{list of files.}
  \item{columns}{what columns to use for signal and background intensities.}
  \item{phenoData}{phenoData object with sample information.}
  \item{\dots}{further arguments to be passed.}
}
\value{readCodelink returns an object of class "Codelink", whereas readCodelinkSet() returns a CodelinkSet object.}
\author{Diego Diez}
\examples{
\dontrun{
    # reading default extension (TXT).
    data <- readCodelink()

    # specify a different one.
    files <- list.files(pattern = "txt")
    data <- readCodelink(files = files)

    # example.
    data(codelink.example)
    
    # read a targets file with phenotypic information
    # including sample names, file names and sample groups.
    pdata = read.AnnotationDataFrame("targets.txt")
    codset = readCodelinkSet(pdata$FileNames, phenoData=pdata)
    
    # example
    data(codelink.exprset)
}
}
\keyword{documentation}
\keyword{utilities}
