\name{ContIC}
\alias{ContIC}

\title{Generating function for ContIC-class}
\description{
  Generates an object of class \code{"ContIC"};
  i.e., an influence curves \eqn{\eta}{eta} of the form
  \deqn{\eta = (A\Lambda - a)\min(1,b/|A\Lambda - a|)}{eta = (A Lambda - a)min(1, b/|A Lambda - a|)}
  with clipping bound \eqn{b}, centering constant \eqn{a} and
  standardizing matrix \eqn{A}. \eqn{\Lambda}{Lambda} stands for
  the L2 derivative of the corresponding L2 differentiable 
  parametric family which can be created via \code{CallL2Fam}.
}
\usage{
ContIC(name, CallL2Fam = call("L2ParamFamily"), 
       Curve = EuclRandVarList(RealRandVariable(Map = c(function(x){x}), 
                                                Domain = Reals())), 
       Risks, Infos, clip = Inf, cent = 0, stand = as.matrix(1), 
       lowerCase = NULL, neighborRadius = 0, w = new("HampelWeight"),
       normtype = NormType(), biastype = symmetricBias(),
       modifyIC = NULL)
}
\arguments{
  \item{name}{ object of class \code{"character"}. }
  \item{CallL2Fam}{ object of class \code{"call"}:
      creates an object of the underlying L2-differentiable 
      parametric family. }
  \item{Curve}{ object of class \code{"EuclRandVarList"} }
  \item{Risks}{ object of class \code{"list"}: 
      list of risks; cf. \code{\link[distrMod]{RiskType-class}}. }
  \item{Infos}{ matrix of characters with two columns 
    named \code{method} and \code{message}: additional informations. }
  \item{clip}{ positive real: clipping bound. }
  \item{cent}{ real: centering constant }
  \item{stand}{ matrix: standardizing matrix }
  \item{w}{ HampelWeight: weight object}
  \item{lowerCase}{ optional constant for lower case solution. }
  \item{neighborRadius}{ radius of the corresponding (unconditional) 
    contamination neighborhood. }
  \item{biastype}{ BiasType: type of the bias}
  \item{normtype}{ NormType: type of the norm}
  \item{modifyIC}{ object of class \code{"OptionalFunction"}:
    function of four arguments: (1) \code{L2Fam} an L2 parametric family
    (2) \code{IC} an optional influence curve, (3) \code{withMakeIC}
    a logical argument whether to enforce the IC side conditions
    by \code{makeIC}, and (4) \code{\dots} for arguments to be passed to
    calls to \code{E} in \code{makeIC}. Returns an object of
    class \code{"IC"}. This function is mainly used for internal
    computations! }
}
%\details{}
\value{Object of class \code{"ContIC"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{IC-class}}, \code{\link{ContIC}} , \code{\link{HampIC-class}}}
\examples{
IC1 <- ContIC()
plot(IC1)
}
\concept{influence curve}
\keyword{robust}
