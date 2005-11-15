\name{fc.Codelink}
\alias{fc.Codelink}
\title{Fold change cutoff calculation}
\description{
	Takes a Codelink object and calculate fold changes between different
	samples. Then look for genes that pass a cutoff.
}
\usage{
	fc.Codelink(object, cond1=NULL, cond2=NULL, fc=1.0)
}
\arguments{
	\item{object}{an object of class "Codelink".}
	\item{cond1}{number of the first condition to use in fold change calculation}
	\item{cond2}{number of the second condition to use in fold change calculation}
	\item{fc}{value of the fold change cutoff}
}
\value{A logical vector indicatin which genes pass the cutoff}
\details{
}
\note{
}
\author{Diego Díez Ruiz}
\seealso{
}
\examples{
}
\keyword{documentation}
\keyword{utilities}
