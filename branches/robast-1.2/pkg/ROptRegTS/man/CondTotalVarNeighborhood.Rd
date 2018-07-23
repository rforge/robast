\name{CondTotalVarNeighborhood}
\alias{CondTotalVarNeighborhood}

\title{Generating function for CondContNeighborhood-class}
\description{
  Generates an object of class \code{"CondTotalVarNeighborhood"}.
}
\usage{
CondTotalVarNeighborhood(radius = 0, radiusCurve = function(x){1})
}
\arguments{
  \item{radius}{ non-negative real: neighborhood radius. }
  \item{radiusCurve}{ real-valued, non-negative function. }
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
\seealso{\code{\link{CondTotalVarNeighborhood-class}}}
\examples{
CondTotalVarNeighborhood()

## The function is currently defined as
function(radius = 0, radiusCurve = function(x){1}){ 
    new("CondTotalVarNeighborhood", radius = radius, radiusCurve = radiusCurve) 
}
}
\concept{conditional total variation neighborhood}
\concept{conditional neighborhood}
\concept{neighborhood}
\keyword{models}
\keyword{robust}
