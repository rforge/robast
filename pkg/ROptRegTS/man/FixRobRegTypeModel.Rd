\name{FixRobRegTypeModel}
\alias{FixRobRegTypeModel}

\title{Generating function for FixRobRegTypeModel-class}
\description{
  Generates an object of class \code{"FixRobRegTypeModel"}.
}
\usage{
FixRobRegTypeModel(center = RegTypeFamily(), neighbor = ContNeighborhood())
}
\arguments{
  \item{center}{ object of class \code{"RegTypeFamily"} }
  \item{neighbor}{ object of class \code{"Neighborhood"} }
}
%\details{}
\value{Object of class \code{"FixRobRegTypeModel"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{FixRobRegTypeModel-class}}}
\examples{
FixRobRegTypeModel()

## The function is currently defined as
function(center = RegTypeFamily(), neighbor = ContNeighborhood()){
    new("FixRobRegTypeModel", center = center, neighbor = neighbor)
}
}
\concept{robust model}
\concept{fixed neighborhood}
\keyword{models}
