\name{checkIC-methods}
\alias{makeIC}
\alias{checkIC}
\alias{makeIC,ContIC,L2ParamFamily-method}
\alias{checkIC,ContIC,L2ParamFamily-method}

\title{Methods for Checking and Making ICs}
\description{
  Particular methods for checking centering and Fisher
  consistency of ICs, resp. making an IC out of an
  IC possibly violating the conditions so far.
}
\usage{
\S4method{checkIC}{ContIC,L2ParamFamily}(IC, L2Fam, out = TRUE,...)
\S4method{makeIC}{ContIC,L2ParamFamily}(IC, L2Fam, ...)
}
\arguments{
  \item{IC}{ object of class \code{"IC"} }
  \item{L2Fam}{ L2-differentiable family of probability measures. }
  \item{out}{ logical: Should the values of the checks be printed out?}
  \item{\dots}{ additional parameters to be passed on to expectation
  \code{E}. }
}
\details{ 
  In \code{checkIC}, the precisions of the centering and the Fisher consistency
  are computed. \code{makeIC} affinely transforms a given IC (not necessarily
  satisfying the centering and Fisher consistency condition so far) such that
  after this transformation it becomes an IC (satisfying the conditions).
  Here particular methods for ICs of class \code{ContIC} are provided using
  the particular structure of this class which allows for speed up in certain cases.
}
\value{The maximum deviation from the IC properties is returned.}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\link{IC-class}}}
\examples{
IC1 <- new("IC")
checkIC(IC1)
}
\concept{influence curve}
\keyword{robust}
