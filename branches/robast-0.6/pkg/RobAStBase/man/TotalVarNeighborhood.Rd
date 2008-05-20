\name{TotalVarNeighborhood}
\alias{TotalVarNeighborhood}

\title{Generating function for TotalVarNeighborhood-class}
\description{
  Generates an object of class \code{"TotalVarNeighborhood"}.
}
\usage{TotalVarNeighborhood(radius = 0)}
\arguments{
  \item{radius}{ non-negative real: neighborhood radius. }
}
%\details{}
\value{Object of class \code{"ContNeighborhood"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{TotalVarNeighborhood-class}}}
\examples{
TotalVarNeighborhood()

## The function is currently defined as
function(radius = 0){ 
    new("TotalVarNeighborhood", radius = radius) 
}
}
\concept{contamination neighborhood}
\concept{neighborhood}
\keyword{models}
