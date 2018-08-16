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
\S4method{checkIC}{ContIC,L2ParamFamily}(IC, L2Fam, out = TRUE,
              forceContICMethod = FALSE, ..., diagnostic = FALSE)
\S4method{makeIC}{ContIC,L2ParamFamily}(IC, L2Fam,
              forceContICMethod = FALSE, ..., diagnostic = FALSE)
}
\arguments{
  \item{IC}{ object of class \code{"IC"} }
  \item{L2Fam}{ L2-differentiable family of probability measures. }
  \item{out}{ logical: Should the values of the checks be printed out?}
  \item{forceContICMethod}{ logical: Should we force to use the method for
  signature \code{ContIC,L2ParamFamily}
  in any case (even if it is not indicated by symmetry arguments)?
  Otherwise it uses internal method \code{.getComp} to compute the number
  of integrals to be computed, taking care of symmetries as indicated through
  the symmetry slots of the model \code{L2Fam}. Only if this
  number is smaller than the  number of integrals to be computed in the range
  of the pIC the present method is used, otherwise it switches back to the
  \code{IC,L2ParamFamily}  method.
  -- The \code{ContIC,L2ParamFamily} up to skipped entries due to further
  symmetry arguments is $\code{(k+1)k/2+k+1=(k+1)(k+2)/2} for \code{k} the
  length of the unknown parameter / length of slot \code{L2deriv} of \code{L2Fam},
  while the number of integrals on the pIC scale underlying the more general
  method for   signature \code{ContIC,L2ParamFamily} is \code{p (k+1)} where
  \code{p} is the length of the pIC / the length of the parameter of interest
  as indicated in the number of rows in the \code{trafo} slot of the underlying
  slot \code{param} of \code{L2Fam}.}
  \item{\dots}{ additional parameters to be passed on to expectation
  \code{E}. }
  \item{diagnostic}{ logical; if \code{TRUE} (and in case \code{checkIC} if
      argument \code{out==TRUE}), diagnostic information on the integration
      is printed and returned as attribute \code{diagnostic} of the return value. }

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
