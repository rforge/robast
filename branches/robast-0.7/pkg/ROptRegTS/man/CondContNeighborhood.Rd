\name{CondContNeighborhood}
\alias{CondContNeighborhood}

\title{Generating function for CondContNeighborhood-class}
\description{
  Generates an object of class \code{"CondContNeighborhood"}.
}
\usage{
CondContNeighborhood(radius = 0, radiusCurve = function(x){1})
}
\arguments{
  \item{radius}{ non-negative real: neighborhood radius. }
  \item{radiusCurve}{ real-valued, non-negative function. }
}
%\details{}
\value{Object of class \code{"CondContNeighborhood"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondContNeighborhood-class}}}
\examples{
CondContNeighborhood()

## The function is currently defined as
function(radius = 0, radiusCurve = function(x){1}){ 
    new("CondContNeighborhood", radius = radius, radiusCurve = radiusCurve) 
}
}
\concept{conditional contamination neighborhood}
\concept{contamination neighborhood}
\concept{neighborhood}
\keyword{models}
