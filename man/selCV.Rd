\name{selCV}
\alias{selCV}
\title{Select based on CV cutoff.}
\description{
	Takes a Codelink object and select genes based on CV cutoff.
}
\usage{
	selCV(object, cutoff)
}
\arguments{
	\item{object}{an object of class "Codelink".}
	\item{cutoff}{cutoff normaly calculated with cutCV()}
}
\value{
	A logical vector.
}
\author{Diego Diez}
\keyword{documentation}
\keyword{utilities}
