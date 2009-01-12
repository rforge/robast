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
}
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[base]{options}}, \code{\link[base]{getOption}}}
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
