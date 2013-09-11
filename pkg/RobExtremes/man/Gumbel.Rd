\name{Gumbel}
\alias{Gumbel}

\title{Generating function for Gumbel-class}
\description{
  Generates an object of class \code{"Gumbel"}.
}
\usage{Gumbel(loc = 0, scale = 1)}
\arguments{
  \item{loc}{ real number: location parameter of 
    the Gumbel distribution. }
  \item{scale}{ positive real number: scale parameter
    of the Gumbel distribution }
}
%\details{}
\value{Object of class \code{"Gumbel"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\note{The class \code{"Gumbel"} is based on the code provided 
  by the package \pkg{evd}.}
\seealso{\code{\link{Gumbel-class}}, \code{\link[evd:gumbel]{rgumbel}}}
\examples{
(G1 <- Gumbel(loc = 1, scale = 2))
plot(G1)
loc(G1)
scale(G1)
loc(G1) <- -1
scale(G1) <- 2
plot(G1)

E(Gumbel()) # Euler's constant
E(G1, function(x){x^2})

## The function is currently defined as
function(loc = 0, scale = 1){ 
  new("Gumbel", loc = loc, scale = scale)
}
}
\concept{Gumbel}
\keyword{distribution}
\concept{absolutely continuous distribution}
\concept{Gumbel distribution}
\concept{generating function}
