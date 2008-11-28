\name{NormLinRegScaleFamily}
\alias{NormLinRegScaleFamily}

\title{ Generating Function for Linear Regression Family with Unknown Scale }
\description{
  Generates an object of class \code{"L2RegTypeFamily"}
  which represents a linear regression family with standard
  normal distributed errors and random regressor where the
  scale of the error distribution is unknown.
}
\usage{
NormLinRegScaleFamily(theta, scale = 1, trafo, RegDistr = Norm(), 
                      RegSymm, Reg2Mom, nuisance = FALSE)
}
\arguments{
  \item{theta}{ linear regression parameter }
  \item{scale}{ scale parameter for error distribution }
  \item{trafo}{ matrix: transformation of the parameter }
  \item{RegDistr}{ regressor distribution }
  \item{RegSymm}{ symmetry of the regressor distribution }
  \item{Reg2Mom}{ second moment matrix of regressor }
  \item{nuisance}{ logical: is scale nuisance parameter }
}
\details{
  In case \code{theta} is missing, it is set to 0.
  If \code{Reg2Mom} is missing, it is computed via \code{E}.
}
\value{Object of class \code{"L2RegTypeFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the
  Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2RegTypeFamily-class}}}
\examples{
(LM1 <- NormLinRegScaleFamily(Reg2Mom = matrix(1)))
Map(L2deriv(LM1)[[1]])
FisherInfo(LM1)
checkL2deriv(LM1)
}
\concept{linear regression}
\concept{linear model}
\keyword{models}
\keyword{robust}
