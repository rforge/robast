\name{InfRobModel}
\alias{InfRobModel}

\title{Generating function for InfRobModel-class}
\description{
  Generates an object of class \code{"InfRobModel"}.
}
\usage{
InfRobModel(center = L2ParamFamily(), neighbor = ContNeighborhood())
}
\arguments{
  \item{center}{ object of class \code{"ProbFamily"} }
  \item{neighbor}{ object of class \code{"UncondNeighborhood"} }
}
%\details{}
\value{Object of class \code{"FixRobModel"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{RobModel-class}}, \code{\link{FixRobModel-class}}}
\examples{
(M1 <- InfRobModel())

## The function is currently defined as
function(center = L2ParamFamily(), neighbor = ContNeighborhood()){
    new("InfRobModel", center = center, neighbor = neighbor)
}
}
\concept{robust model}
\concept{infinitesimal neighborhood}
\concept{shrinking neighborhood}
\keyword{models}
