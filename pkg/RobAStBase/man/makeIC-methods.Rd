\name{makeIC}
%\docType{methods}
\alias{makeIC}
\alias{makeIC-methods}
\alias{makeIC,IC,missing-method}
\alias{makeIC,IC,L2ParamFamily-method}
\alias{makeIC,list,L2ParamFamily-method}
\alias{makeIC,function,L2ParamFamily-method}

\title{Generic Function for making ICs consistent at a possibly different model}
\description{
  Generic function for providing centering and Fisher consistency of ICs.
}
\usage{
makeIC(IC, L2Fam, ...)
%\S4method{makeIC}{IC,missing}(IC, ...)
\S4method{makeIC}{IC,L2ParamFamily}(IC, L2Fam, ..., diagnostic = FALSE)
\S4method{makeIC}{list,L2ParamFamily}(IC, L2Fam, forceIC = TRUE, name, Risks,
                  Infos, modifyIC = NULL, ..., diagnostic = FALSE)
\S4method{makeIC}{function,L2ParamFamily}(IC, L2Fam, forceIC = TRUE, name, 
                  Risks, Infos, modifyIC = NULL, ..., diagnostic = FALSE)
}
\arguments{
  \item{IC}{ object of class \code{"IC"} for signature \code{IC="IC"}, respectively
   a list of functions in one argument for signature \code{IC="list"}, respectively
   a function in one argument for signature \code{IC="function"}.}
  \item{L2Fam}{ L2-differentiable family of probability measures; may be missing,
                in which case it is replaced by the family in slot \code{CallL2Fam}
                of \code{IC}. }
  \item{forceIC}{ logical; shall centeredness and Fisher consistency be enforced
                  applying an affine linear transformation?}
  \item{name}{ Object of class \code{"character"}; the name of the IC }
  \item{Risks}{ object of class \code{"list"}: 
      list of risks; cf. \code{\link[distrMod]{RiskType-class}}. }
  \item{Infos}{ matrix of characters with two columns 
    named \code{method} and \code{message}: additional informations. }
  \item{modifyIC}{ object of class \code{"OptionalFunction"}:
    function of four arguments: (1) \code{L2Fam} an L2 parametric family
    (2) \code{IC} an optional influence curve, (3) \code{withMakeIC}
    a logical argument whether to enforce the IC side conditions
    by \code{makeIC}, and (4) \code{\dots} for arguments to be passed to
    calls to \code{E} in \code{makeIC}. Returns an object of
    class \code{"IC"}. This function is mainly used for internal
    computations! }
  \item{\dots}{ additional parameters to be passed to expectation \code{E} }
  \item{diagnostic}{ logical; if \code{TRUE},
    diagnostic information on the integration is printed and returned
    as attribute \code{diagnostic} of the return value. }
}
\value{An IC of class \code{"IC"} at the model.}
\details{
  Argument \code{IC} is transformed affinely such that the transformed IC
  satisfies the defining side conditions of an IC, i.e., centeredness and
  Fisher consistency:

  \deqn{\mathop{\mbox{\boldmath$E$}}[{\rm IC}]=0}{E[IC]=0}
  \deqn{\mathop{\mbox{\boldmath$E$}}[{\rm IC}\,\Lambda^\tau]= D}{E[IC Lambda'] = D}

  where \eqn{\Lambda}{Lambda} is the L2 derivative of the model and D is
  the Jacobian of transformation \code{trafo}.

  Diagnostics on the involved integrations are available if argument
  \code{diagnostic} is \code{TRUE}. Then there is attribute \code{diagnostic}
  attached to the return value, which may be inspected
  and accessed through \code{\link[distrEx:distrExIntegrate]{showDiagnostic}} and
  \code{\link[distrEx:distrExIntegrate]{getDiagnostic}}.
}

\section{Methods}{\describe{
\item{makeIC}{\code{signature(IC = "IC", L2Fam = "missing"}: creates 
          an object of class \code{"IC"} at the parametric model of its own
          slot \code{CallL2Fam}; enforces IC conditions 
          centeredness and Fisher consistency, applying an affine linear 
          transformation.}
\item{makeIC}{\code{signature(IC = "IC", L2Fam = "L2ParamFamily"}: creates 
          an object of class \code{"IC"} at the parametric model \code{L2Fam}; 
          enforces IC conditions centeredness and Fisher consistency, 
          applying an affine linear transformation.}
\item{makeIC}{\code{signature(IC = "list", L2Fam = "L2ParamFamily"}: creates 
          an object of class \code{"IC"} out of a list of functions given by argument 
          \code{IC} at the parametric model \code{L2Fam}; 
          enforces IC conditions centeredness and Fisher consistency, 
          applying an affine linear transformation.}
\item{makeIC}{\code{signature(IC = "function", L2Fam = "L2ParamFamily"}: creates 
          an object of class \code{"IC"} out of a function given by argument 
          \code{IC} at the parametric model \code{L2Fam}; 
          enforces IC conditions centeredness and Fisher consistency, 
          applying an affine linear transformation.}
}}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
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

IC4 <- makeIC(sin, B)
checkIC(IC4)

(IC5 <- makeIC(list(function(x)x^3), B, name="a try"))
plot(IC5)
checkIC(IC5)

## don't run to reduce check time on CRAN
\donttest{
N0 <- NormLocationScaleFamily()
IC6 <- makeIC(list(sin,cos),N0)
plot(IC6)
checkIC(IC6)

getRiskIC(IC6,risk=trAsCov())$trAsCov$value
getRiskIC(IC6,risk=asBias(),neighbor=ContNeighborhood())$asBias$value
}

}
\concept{influence curve}
\keyword{robust}
