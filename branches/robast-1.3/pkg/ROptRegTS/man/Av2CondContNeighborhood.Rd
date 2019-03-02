\name{Av2CondContNeighborhood}
\alias{Av2CondContNeighborhood}

\title{Generating function for Av2CondContNeighborhood-class}
\description{
  Generates an object of class \code{"Av2CondContNeighborhood"}.
}
\usage{
Av2CondContNeighborhood(radius = 0, radiusCurve = function(x){1})
}
\arguments{
  \item{radius}{ non-negative real: neighborhood radius. }
  \item{radiusCurve}{ real-valued, non-negative function
    with L2 norm \code{<= 1}. }
}
%\details{}
\value{Object of class \code{"Av1CondContNeighborhood"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{Av2CondContNeighborhood-class}}}
\examples{
Av2CondContNeighborhood()

## The function is currently defined as
function(radius = 0, radiusCurve = function(x){1}){ 
    new("Av2CondContNeighborhood", radius = radius, radiusCurve = radiusCurve) 
}
}
\concept{average conditional neighborhood}
\concept{conditional neighborhood}
\concept{neighborhood}
\keyword{models}
\keyword{robust}
