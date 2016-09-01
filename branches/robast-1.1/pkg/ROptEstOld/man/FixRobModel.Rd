\name{FixRobModel}
\alias{FixRobModel}

\title{Generating function for FixRobModel-class}
\description{
  Generates an object of class \code{"FixRobModel"}.
}
\usage{
FixRobModel(center = ParamFamily(), neighbor = ContNeighborhood())
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
\seealso{\code{\link{FixRobModel-class}}}
\examples{
(M1 <- FixRobModel())

## The function is currently defined as
function(center = ParamFamily(), neighbor = ContNeighborhood()){
    new("FixRobModel", center = center, neighbor = neighbor)
}
}
\concept{robust model}
\concept{fixed neighborhood}
\keyword{models}
