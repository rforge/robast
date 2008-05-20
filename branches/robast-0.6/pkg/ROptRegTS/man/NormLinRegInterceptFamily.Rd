\name{NormLinRegInterceptFamily}
\alias{NormLinRegInterceptFamily}

\title{ Generating Function for Linear Regression Family with Unknown Intercept }
\description{
  Generates an object of class \code{"L2RegTypeFamily"}
  which represents a linear regression family with standard
  normal distributed errors and random regressor where the
  intercept is unknown.
}
\usage{
NormLinRegInterceptFamily(theta, intercept = 0, trafo, RegDistr = Norm(), 
                          RegSymm, Reg2Mom, nuisance = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{ linear regression parameter }
  \item{intercept}{ intercept parameter }
  \item{trafo}{ matrix: transformation of the parameter }
  \item{RegDistr}{ regressor distribution }
  \item{RegSymm}{ symmetry of the regressor distribution }
  \item{Reg2Mom}{ second moment matrix of regressor }
  \item{nuisance}{ logical: is intercept nuisance parameter }
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
(LM1 <- NormLinRegInterceptFamily(Reg2Mom = matrix(1)))
Map(L2deriv(LM1)[[1]])
FisherInfo(LM1)
checkL2deriv(LM1)
}
\concept{linear regression}
\concept{linear model}
\keyword{models}
