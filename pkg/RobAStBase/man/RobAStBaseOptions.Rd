\name{RobAStBaseOptions}
\alias{RobAStBaseOptions}
\alias{getRobAStBaseOption}
\alias{kStepUseLast}

\title{Function to change the global variables of the package `RobAStBase' }
\description{With \code{RobAStBaseOptions} you can inspect and change 
  the global variables of the package \pkg{RobAStBase}. }
\usage{
RobAStBaseOptions(...)
getRobAStBaseOption(x)
}
\arguments{
  \item{\dots}{ any options can be defined, using name = value or by passing a list of such tagged values. }
  \item{x}{ a character string holding an option name.}
}
%\details{}
\value{
  \code{RobAStBaseOptions()} returns a list of the global variables.\cr
  \code{RobAStBaseOptions(x)} returns the global variable \var{x}.\cr
  \code{getRobAStBaseOption(x)} returns the global variable \var{x}.\cr
  \code{RobAStBaseOptions(x=y)} sets the value of the global variable \var{x} to \var{y}.
}
\section{Global Options}{
\describe{
  \item{kStepUseLast:}{ The default value of argument \code{kStepUseLast} is
    \code{FALSE}. Explicitly setting \code{kStepUseLast} to \code{TRUE} should
    be done with care as in this situation the influence curve in case of
    \code{oneStepEstimator} and \code{kStepEstimator} is re-computed using 
    the value of the one- resp. k-step estimate which may take quite a long 
    time depending on the model. }
  \item{withUpdateInKer:}{if there is a non-trivial trafo in the model with matrix \eqn{D}, shall
  the parameter be updated on \eqn{{\rm ker}(D)}{ker(D)}? defaults to \code{FALSE}.}
  \item{IC.UpdateInKer:}{if there is a non-trivial trafo in the model with matrix \eqn{D},
  the IC to be used for this; if \code{NULL} the result of \code{getboundedIC(L2Fam,D)} is taken;
  this IC will then be projected onto \eqn{{\rm ker}(D)}{ker(D)}; defaults to \code{NULL}.}
  \item{all.verbose:}{argument \code{verbose} passed on by default to many calls of
  \code{optIC}, \code{radiusminimaxIC}, \code{getinfRobIC} etc.; well suited
  for testing purposes.}
  \item{withPICList:}{logical: shall slot \code{pICList} of return value of \code{kStepEstimator}
  be filled?}
  \item{withICList:}{logical: shall slot \code{ICList} of return value of \code{kStepEstimator}
  be filled?}
}
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[base]{options}}, \code{\link[base:options]{getOption}}}
\examples{
RobAStBaseOptions()
RobAStBaseOptions("kStepUseLast")
RobAStBaseOptions("kStepUseLast" = TRUE)
# or
RobAStBaseOptions(kStepUseLast = 1e-6)
getRobAStBaseOption("kStepUseLast")
}
\keyword{misc}
\keyword{robust}
\concept{global options}
\concept{options}
