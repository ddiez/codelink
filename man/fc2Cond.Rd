\name{fc2Cond}
\alias{fc2Cond}
\title{Select probes based on fold change calculation}
\description{
  Takes a Codelink object and calculate fold changes (M) between two conditions 
(samples). Then select genes based on those who pass the passed cutoff.
}
\usage{
  fc2Cond(object, cond1=NULL, cond2=NULL, fc=1.0, verbose=FALSE)
}
\arguments{
  \item{object}{an object of class "Codelink".}
  \item{cond1}{numeric or character; First condition to compute M.}
  \item{cond2}{numeric or character; Second condition to compute M.}
  \item{fc}{value of the fold change cutoff}
  \item{verbose}{logical; if some information is printed on the console.}
}
\value{A logical vector indicating which genes pass the cutoff}
\details{
  Conditions can be passed as characters or as numeric index from the sample 
slot. The intensities are internally transformed to log2 if needed. The M value 
is computed as:

   M = cond1 - cond2
}
\author{Diego Diez}
\keyword{documentation}
\keyword{utilities}
