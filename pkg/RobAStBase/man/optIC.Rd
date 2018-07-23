\name{optIC}
\alias{optIC}
\alias{optIC-methods}
\alias{optIC,L2ParamFamily,asCov-method}

\title{Generic function for the computation of optimally robust ICs}
\description{
  Generic function for the computation of optimally robust ICs.
}
\usage{
optIC(model, risk, ...)

\S4method{optIC}{L2ParamFamily,asCov}(model, risk, withMakeIC = FALSE)
}
\arguments{
  \item{model}{ probability model. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{\dots}{ additional parameters. }
  \item{withMakeIC}{logical; if \code{TRUE} the [p]IC is passed through
        \code{makeIC} before return.}
}
\details{ 
  The classical optimal IC which ist optimal in sense of the Cramer-Rao bound
  is computed.
}
\value{Some optimally robust IC is computed.}
\section{Methods}{
\describe{
  \item{model = "L2ParamFamily", risk = "asCov"}{ computes
    classical optimal influence curve for L2 differentiable 
    parametric families.}
}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{InfluenceCurve-class}}, \code{\link[distrMod]{RiskType-class}}}
\examples{
B <- BinomFamily(size = 25, prob = 0.25) 

## classical optimal IC
IC0 <- optIC(model = B, risk = asCov())
plot(IC0) # plot IC
checkIC(IC0, B)
}
\concept{robust influence curve}
\concept{influence curve}
\keyword{robust}
