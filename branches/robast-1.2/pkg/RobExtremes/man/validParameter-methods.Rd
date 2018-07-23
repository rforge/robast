\name{validParameter-methods}
\docType{methods}
\alias{validParameter-methods}
\alias{validParameter}
\alias{validParameter,GParetoFamily-method}
\alias{validParameter,WeibullFamily-method}
\alias{validParameter,GEVFamily-method}
\alias{validParameter,ParetoFamily-method}
\alias{validParameter,GEVFamilyMuUnknown-method}


\title{ Methods for function validParameter in Package `RobExtremes' }
\description{
  Methods for function \code{validParameter} in package \pkg{RobExtremes}
  to check whether a new parameter (e.g. "proposed" by an optimization)
  is valid.
}

\usage{
validParameter(object, \dots)
\S4method{validParameter}{GParetoFamily}(object, param, tol=.Machine$double.eps)
\S4method{validParameter}{WeibullFamily}(object, param, tol=.Machine$double.eps)
\S4method{validParameter}{GEVFamily}(object, param, tol=.Machine$double.eps)
\S4method{validParameter}{ParetoFamily}(object, param, tol=.Machine$double.eps)
\S4method{validParameter}{GEVFamilyMuUnknown}(object, param,
           tol=.Machine$double.eps)
}

\details{
method for signature
\describe{
  \item{\code{GParetoFamily}}{checks if both parameters are finite by \code{is.finite},
  if their length is 1 or 2 (e.g. if one features as nuisance parameter), and if
  both are strictly larger than 0 (upto argument \code{tol})}
  \item{\code{WeibullFamily}}{checks if both parameters are finite by \code{is.finite},
  if their length is 1 or 2 (e.g. if one features as nuisance parameter), and if
  both are strictly larger than 0 (upto argument \code{tol})}
  \item{\code{GEVFamily}}{checks if both parameters are finite by \code{is.finite},
  if their length is 1 or 2 (e.g. if one features as nuisance parameter), and if
  both are strictly larger than 0 (upto argument \code{tol})}
  \item{\code{GParetoFamily}}{checks if both parameters are finite by \code{is.finite},
  if their length is 1 or 2 (e.g. if one features as nuisance parameter), and if
  both are strictly larger than 0 (upto argument \code{tol})}
  \item{\code{GEVFamilyMuUnknown}}{checks if all parameters are finite by \code{is.finite},
  if their length is in 1,2,3 (e.g. if one features as nuisance parameter), and scale
  and shape both are strictly larger than 0 (upto argument \code{tol})}
}
}
\arguments{
  \item{object}{an object of class \code{ParamFamily}}
  \item{param}{either a numeric vector or an object of class 
  \code{ParamFamParameter}}
  \item{tol}{accuracy upto which the conditions have to be fulfilled}
  \item{\dots}{additional argument(s) for methods.}
}
\value{\code{logical} of length 1 --- valid or not}
\examples{
 G <- GParetoFamily()
 validParameter(G, c(scale=0.1, shape=2))
 validParameter(G, c(scale=-0.1, shape=-2))
}
\keyword{models}
