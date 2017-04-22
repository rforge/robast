\name{GPareto}
\alias{GPareto}

\title{Generating function for GPareto-class}
\description{
  Generates an object of class \code{"GPareto"}.
}
\usage{GPareto(loc = 0, scale = 1, shape = 0, location = loc)}
\arguments{
  \item{loc}{ real number: location parameter of 
    the GPareto distribution. }
   \item{scale}{ positive real number: scale parameter
    of the GPareto distribution }
   \item{shape}{ non-negative real number: shape parameter of 
    the GPareto distribution. }
  \item{location}{ alternative argument name for argument 'loc' --- to 
    support argument names of package \pkg{VGAM}. }
 
}
%\details{}
\value{Object of class \code{"GPareto"}}
%\references{}
\author{Nataliya Horbenko \email{nhorbenko@gmail.com}}
\note{The class \code{"GPareto"} is based on the code provided 
  by the package \pkg{evd} by  Alec Stephenson.}
\seealso{\code{\link{GPareto-class}}, \code{\link[evd:gpd]{dgpd}}}
\examples{
(P1 <- GPareto(loc = 1, scale = 1, shape = -0.5))
plot(P1)

E(GPareto()) 
E(P1)
E(P1, function(x){x^2})
var(P1)
sd(P1)
median(P1)
IQR(P1)
mad(P1)

}

\concept{GPareto}
\keyword{distribution}
\concept{absolutely continuous distribution}
\concept{GPareto distribution}
\concept{generating function}
