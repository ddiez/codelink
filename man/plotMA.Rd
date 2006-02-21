\name{plotMA}
\alias{plotMA}
\title{MA plot}
\description{
	Takes a Codelink object and plot M vs A.
}
\usage{
  plotMA(object, array1=1, array2=2, cutoff=NULL, label="type",
	high.list=NULL, high.col="blue", high.pch="*",
	snr.cutoff=1, legend.x="bottomright", pch=".", subset=NULL, 
	title=NULL, xlim=NULL, ylim=NULL)
}
\arguments{
  \item{object}{an object of class "Codelink".}
  \item{array1}{first array to be used.}
  \item{array2}{second array to be used.}
  \item{cutoff}{cutoff to be used as fold change markeer.}
  \item{label}{type of labeling used in legend.}
  \item{high.list}{list of genes highlighted.}
  \item{high.col}{color used for high genes.}
  \item{high.pch}{pch used for high genes.}
  \item{snr.cutoff}{SNR cutoff used for label spots.}
  \item{legend.x}{relative position of the legend.}
  \item{pch}{pch style used to main spots.}
  \item{subset}{subset of spots used to plot based on 'type' slot.}
  \item{title}{title of the plot.}
  \item{xlim}{range for the X axis.}
  \item{ylim}{range for the Y axis.}
}
\value{
}
\details{
	The plot is for M=(array2-array1) vs A=(array2+array1)/2.
}
\note{
}
\author{Diego Díez Ruiz}
\seealso{
}
\examples{
#    data(codelink.example)
#    plotMA(codelink.example)
}
\keyword{documentation}
\keyword{utilities}
