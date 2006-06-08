\name{na2false}
\alias{na2false}
\title{Set NAs to FALSE}
\description{
	Takes a logical vector as input and set all NAs to FALSE. This may happens when comparison is done on NA values.
}
\usage{
	na2false(x)
}
\arguments{
  \item{x}{a logical vector.}
}
\value{A logical vector without NAs.}
\details{
}
\note{
}
\author{Diego Diez}
\seealso{
}
\examples{
 #   a <- c(1, 2, 3, NA, 5)
 #   b <- c(5, 4, NA, 2, 1)
 #   sel <- a > b
 #   sel <- na2false(sel)
}
\keyword{documentation}
\keyword{utilities}
