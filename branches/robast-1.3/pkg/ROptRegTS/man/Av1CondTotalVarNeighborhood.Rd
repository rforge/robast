\name{Av1CondTotalVarNeighborhood}
\alias{Av1CondTotalVarNeighborhood}

\title{Generating function for Av1CondTotalVarNeighborhood-class}
\description{
  Generates an object of class \code{"Av1CondTotalVarNeighborhood"}.
}
\usage{
Av1CondTotalVarNeighborhood(radius = 0, radiusCurve = function(x){1})
}
\arguments{
  \item{radius}{ non-negative real: neighborhood radius. }
  \item{radiusCurve}{ real-valued, non-negative function
    with L1 norm \code{<= 1}. }
}
%\details{}
\value{Object of class \code{"Av1CondTotalVarNeighborhood"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{Av1CondTotalVarNeighborhood-class}}}
\examples{
Av1CondTotalVarNeighborhood()

## The function is currently defined as
function(radius = 0, radiusCurve = function(x){1}){ 
    new("Av1CondTotalVarNeighborhood", radius = radius, radiusCurve = radiusCurve) 
}
}
\concept{average conditional neighborhood}
\concept{conditional neighborhood}
\concept{neighborhood}
\keyword{models}
\keyword{robust}
