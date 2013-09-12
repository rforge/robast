\name{infoPlot}
\alias{infoPlot}

\title{Plot absolute and relative information}
\description{
  Plot absolute and relative information of influence curves.
}
\usage{
infoPlot(object)
}
\arguments{
  \item{object}{ object of class \code{"InfluenceCurve"} }
}
\details{
  Absolute information is defined as the square of the length
  of an IC. The relative information is defined as the 
  absolute information of one component with respect to the 
  absolute information of the whole IC; confer Section 8.1 
  of Kohl (2005).
}
%\value{}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link{IC-class}}}
\examples{
N <- NormLocationScaleFamily(mean=0, sd=1) 
IC1 <- optIC(model = N, risk = asCov())
infoPlot(IC1)
}
\concept{absolute information}
\concept{relative information}
\keyword{robust}
