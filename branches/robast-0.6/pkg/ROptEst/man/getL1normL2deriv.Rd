\name{getL1normL2deriv}
\alias{getL1normL2deriv}
\alias{getL1normL2deriv-methods}
\alias{getL1normL2deriv,UnivariateDistribution-method}
\alias{getL1normL2deriv,RealRandVariable-method}

\title{Calculation of L1 norm of L2derivative}
\description{
  Methods to calculate the L1 norm of the L2derivative in a smooth parametric model.
}
\usage{getL1normL2deriv(L2deriv, ...)
\S4method{getL1normL2deriv}{UnivariateDistribution}(L2deriv, 
     cent, ...)

\S4method{getL1normL2deriv}{RealRandVariable}(L2deriv, 
     cent, stand, Distr, normtype, ...)
}
%\details{}
\arguments{
  \item{L2deriv}{L2derivative of the model}
  \item{cent}{centering Lagrange Multiplier}
  \item{stand}{standardizing Lagrange Multiplier}
  \item{Distr}{distribution of the L2derivative}
  \item{normtype}{object of class \code{NormType}; the norm under which we work}
  \item{\dots}{further arguments (not used at the moment)}
}

\value{L1 norm of the L2derivative}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\examples{
##
}
\concept{L1norm}
\keyword{robust}

