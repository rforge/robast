\name{makeIC-methods}
\docType{methods}
\alias{makeIC}
\alias{makeIC-methods}
\alias{makeIC,IC,missing-method}
\alias{makeIC,IC,L2ParamFamily-method}

\title{Generic Function for making ICs consistent at a possibly different model}
\description{
  Generic function for providing centering and Fisher consistency of ICs.
}
\usage{
makeIC(IC, L2Fam, ...)
}
\arguments{
  \item{IC}{ object of class \code{"IC"} }
  \item{L2Fam}{ L2-differentiable family of probability measures; may be missing. }
  \item{\dots}{ additional parameters }
}
\value{An IC at the model.}
\section{Methods}{\describe{
\item{makeIC}{\code{signature(IC = "IC", L2Fam = "missing"}: ...}
\item{makeIC}{\code{signature(IC = "IC", L2Fam = "L2ParamFamily"}: ...}
}}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\link{IC-class}}}
\examples{
## default IC
IC1 <- new("IC")

## L2-differentiable parametric family
B <- BinomFamily(13, 0.3)

## check IC properties
checkIC(IC1, B)

## make IC
IC2 <- makeIC(IC1, B)

## check IC properties
checkIC(IC2)

## slot modifyIC is filled in case of IC2
IC3 <- modifyIC(IC2)(BinomFamily(13, 0.2), IC2)
checkIC(IC3)
## identical to
checkIC(IC3, BinomFamily(13, 0.2))
}
\concept{influence curve}
\keyword{robust}
