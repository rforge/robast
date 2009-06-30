\name{InfRobRegTypeModel}
\alias{InfRobRegTypeModel}

\title{Generating function for InfRobRegTypeModel-class}
\description{
  Generates an object of class \code{"InfRobRegTypeModel"}.
}
\usage{
InfRobRegTypeModel(center = L2RegTypeFamily(), neighbor = ContNeighborhood())
}
\arguments{
  \item{center}{ object of class \code{"L2RegTypeFamily"} }
  \item{neighbor}{ object of class \code{"Neighborhood"} }
}
%\details{}
\value{Object of class \code{"InfRobRegTypeModel"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{InfRobRegTypeModel-class}}}
\examples{
InfRobRegTypeModel()

## The function is currently defined as
function(center = L2RegTypeFamily(), neighbor = ContNeighborhood()) {
    new("InfRobRegTypeModel", center = center, neighbor = neighbor)
}
}
\concept{robust model}
\concept{infinitesimal neighborhood}
\concept{shrinking neighborhood}
\keyword{models}
\keyword{robust}
