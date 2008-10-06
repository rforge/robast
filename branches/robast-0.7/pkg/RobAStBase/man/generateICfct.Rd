\name{generateIC.fct-methods}
\docType{methods}
\alias{generateIC.fct}
\alias{generateIC.fct-methods}
\alias{generateIC.fct,UncondNeighborhood,L2ParamFamily-method}

\title{Generic Function for making ICs consistent at a possibly different model}
\description{
  Generic function for providing centering and Fisher consistency of ICs.
}
\usage{
generateIC.fct(neighbor, L2Fam, ...)
}
\arguments{
  \item{neighbor}{ object of class \code{"UncondNeighborhood"} }
  \item{L2Fam}{ L2-differentiable family of probability measures; may be missing. }
  \item{\dots}{ additional parameters }
}
\value{An IC at the model.}
\section{Methods}{\describe{
\item{generateIC.fct}{\code{signature(IC = "UncondNeighborhood", L2Fam = "L2ParamFamily"}: ...}
}}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}.
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\link{IC-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
