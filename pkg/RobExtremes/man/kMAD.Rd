\name{kMAD}
\alias{kMAD}
\alias{kMAD-methods}
\alias{kMAD,UnivariateDistribution,numeric-method}
\alias{kMAD,numeric,numeric-method}

\title{Asymmetric Median of Absolute Deviations for Skewed Distributions}
\description{
  Function for the computation of asymmetric median absolute deviation (kMAD)
  It coincides with ordinary median absolute deviation (MAD) for \eqn{k=1}{k=1}.
}
\usage{
kMAD(x,k,...)
\S4method{kMAD}{numeric,numeric}(x, k = 1, na.rm = TRUE, 
                eps = .Machine$double.eps, ... )
\S4method{kMAD}{UnivariateDistribution,numeric}(x, k = 1, up = NULL, ... )
}
\arguments{
  \item{x}{a numeric vector or a distribution. }
  \item{k}{numeric; tunning parameter for asymmetrical MAD; has to be of length 1 and larger than 1.}
  \item{na.rm}{logical; if \code{TRUE} then \code{NA} values are stripped from \code{x} before computation takes place. }
  \item{eps}{numeric; accuracy up to which to state equality of two numeric values }
  \item{up}{numeric; upper bound for search interval; important in distributions without left/right endpoint.}
  \item{\dots}{additional arguments for other functions; not used so far;}
 }
\details{
  For kMAD (asymmetrial MAD)  is a root of the equation: 
  \deqn{\mathop{\rm kMAD}(F,k) = \inf\{t>0\;\mid \;F(m+kt)-F(m-t)\ge 1/2 \}}{kMAD(F,k) = inf{t>0|F(m+kt)-F(m-t)>=1/2}}, 
  where \code{F} is the cumulative distribution function, \code{m} is the median of \code{F}.
   }
\references{
  Ruckdeschel, P., Horbenko, N. (2010): Robustness Properties for Generalized Pareto Distributions. ITWM Report 182.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}, 
Nataliya Horbenko \email{Nataliya.Horbenko@itwm.fraunhofer.de}}
\seealso{\code{\link{mad}}}
\examples{
x <- rnorm(100)
kMAD(x,k=10)
kMAD(Norm(),k=10)
}
\concept{estimator}
\keyword{scale estimator}
