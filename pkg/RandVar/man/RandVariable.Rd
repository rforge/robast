\name{RandVariable}
\alias{RandVariable}

\title{Generating function for RandVariable-class}
\description{
  Generates an object of class \code{"RandVariable"}.
}
\usage{
RandVariable(Map = list(function(x){}), Domain = NULL, Range = NULL)
}
\arguments{
  \item{Map}{ list of functions forming the map. }
  \item{Domain}{ domain of \code{Map}: object of 
    class \code{"OptionalrSpace"} (default = \code{NULL}). }
  \item{Range}{ range of \code{Map}: object of 
    class \code{"OptionalrSpace"} (default = \code{NULL}). }
}
%\details{}
\value{Object of class \code{"RandVariable"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{RandVariable-class}}}
\examples{
(R1 <- RandVariable())
Map(R1)
Domain(R1)
Range(R1)
Map(R1) <- list(function(x){ceiling(x)}, function(x){floor(x)})
Domain(R1) <- Reals()
Range(R1) <- Naturals()
R1
Map(R1)
length(R1)

R2 <- R1
Domain(R2) <- Naturals()
compatibleDomains(R1, R2)
Domain(R2) <- NULL
compatibleDomains(R1, R2)
Domain(R2) <- EuclideanSpace(dimension = 1)
compatibleDomains(R1, R2)
Domain(R2) <- EuclideanSpace(dimension = 2)
compatibleDomains(R1, R2)

## The function is currently defined as
function(Map = list(function(x){ }), Domain = NULL, Range = NULL) {
    return(new("RandVariable", Map = Map, Domain = Domain, Range = Range))
}
}
\concept{random variable}
\keyword{}
