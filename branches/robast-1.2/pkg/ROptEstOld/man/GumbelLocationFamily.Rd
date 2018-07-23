\name{GumbelLocationFamily}
\alias{GumbelLocationFamily}

\title{Generating function for Gumbel location families}
\description{
  Generates an object of class \code{"L2ParamFamily"} which
  represents a Gumbel location family.
}
\usage{
GumbelLocationFamily(loc = 0, scale = 1, trafo)
}
\arguments{
  \item{loc}{ location parameter }
  \item{scale}{ scale parameter }
  \item{trafo}{ matrix: transformation of the parameter }
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
}
\value{Object of class \code{"L2ParamFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link{Gumbel-class}}}
\examples{
distrExOptions("ElowerTruncQuantile" = 1e-15) # problem with 
                                            # non-finite function value
(G1 <- GumbelLocationFamily())
plot(G1)
Map(L2deriv(G1)[[1]])
checkL2deriv(G1)
distrExOptions("ElowerTruncQuantile" = 0) # default
}
\concept{Gumbel location model}
\concept{location model}
\keyword{models}
