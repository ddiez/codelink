\name{plotMA}
\alias{plotMA}
\title{MA plot}
\description{
	Takes a Codelink object and plot M vs A.
}
\usage{
  plotMA(object, array1 = 1, array2 = NULL, cutoff = c(-1, 1), label = NULL,
    type = NULL, high.list = NULL, high.col = "blue", high.pch = 21,
	high.bg = "cyan",	snr = NULL, snr.cutoff = 1, legend.x = NULL, pch = ".", 
	subset = NULL, title = NULL, xlim = NULL, ylim = NULL)
}
\arguments{
  \item{object}{an object of class "Codelink" or "MAarrayLM".}
  \item{array1}{first array to be used.}
  \item{array2}{second array to be used.}
  \item{cutoff}{cutoff to be used as fold change markeer.}
  \item{label}{type of labeling used in legend.}
  \item{type}{spot type information.}
  \item{high.list}{list of genes highlighted.}
  \item{high.col}{color used for high genes.}
  \item{high.pch}{pch used for high genes.}
  \item{high.bg}{background color used for high genes.}
  \item{snr}{vector with SNR values, usually, taking rowMeans() from a SNR matrix.}
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

  This function has suffered recent re-working, to increase the usability and 
  to clean a little bit the code.

  If array2 is NULL a median array is computed using all available arrays.  
  Then the values of M and A are computed using the following formula:

                             M = array2 - array1

                          A = (array2 + array1) / 2

  If type information is available in the Codelink object, or provided throught
  the 'type' argument, spots are colored based on that. DISCOVERY spots are 
  plotted black with pch = "." whereas the other classes are plotted with 
  different background colors, using gray as border to increase contrasts. For 
  that pch = 21 is used. If snr is specified as label option, the SNR is used 
  to label spots, if available in the Codelink object. In this case, the mean 
  SNR across all arrays is used when array2 = NULL.

  Some parameters may not be working right now, as the new function is using a 
  different method to labels spots.

  The legend is 'automagically' located, but this can be overrided with the 
  legend.x argument.

  In addition, a subset of the spots can be plotted based on type information 
  when available. This allows, for example, to plot only DISCOVERY spots.
}
\note{
}
\author{Diego Diez}
\seealso{
}
\examples{
\dontrun{
    data(codelink.example)
    plotMA(codelink.example)
}
}
\keyword{documentation}
\keyword{utilities}
