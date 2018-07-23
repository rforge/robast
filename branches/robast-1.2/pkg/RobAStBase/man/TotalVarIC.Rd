\name{TotalVarIC}
\alias{TotalVarIC}

\title{Generating function for TotalVarIC-class}
\description{
  Generates an object of class \code{"TotalVarIC"};
  i.e., an influence curves \eqn{\eta}{eta} of the form
  \deqn{\eta = c \vee A\Lambda \wedge d}{eta = max(c, min(A Lambda, d))}
  with lower clipping bound \eqn{c}, upper clipping bound \eqn{d} and
  standardizing matrix \eqn{A}. \eqn{\Lambda}{Lambda} stands for
  the L2 derivative of the corresponding L2 differentiable 
  parametric family which can be created via \code{CallL2Fam}.
}
\usage{
TotalVarIC(name, CallL2Fam = call("L2ParamFamily"), 
           Curve = EuclRandVarList(RealRandVariable(Map = c(function(x) {x}), 
                                                    Domain = Reals())), 
           Risks, Infos, clipLo = -Inf, clipUp = Inf, stand = as.matrix(1), 
           lowerCase = NULL, neighborRadius = 0, w = new("BdStWeight"),
           normtype = NormType(), biastype = symmetricBias(),
           modifyIC = NULL)
}
\arguments{
  \item{name}{ object of class \code{"character"}. }
  \item{CallL2Fam}{ object of class \code{"call"}:
      creates an object of the underlying L2-differentiable 
      parametric family. }
  \item{Curve}{ object of class \code{"EuclRandVarList"}. }
  \item{Risks}{ object of class \code{"list"}: 
      list of risks; cf. \code{\link[distrMod]{RiskType-class}}. }
  \item{Infos}{ matrix of characters with two columns 
    named \code{method} and \code{message}: additional informations. }
  \item{clipLo}{ negative real: lower clipping bound. }
  \item{clipUp}{ positive real: lower clipping bound. }
  \item{stand}{ matrix: standardizing matrix }
  \item{w}{ BdStWeight: weight object}
  \item{lowerCase}{ optional constant for lower case solution. }
  \item{neighborRadius}{ radius of the corresponding (unconditional) 
    contamination neighborhood. }
  \item{biastype}{ BiasType: type of the bias}
  \item{normtype}{ NormType: type of the norm}
  \item{modifyIC}{ object of class \code{"OptionalFunction"}:
    function of two arguments, which are an L2 parametric family 
    and an optional influence curve. Returns an object of 
    class \code{"IC"}. This function is mainly used for internal
    computations! }
}
%\details{}
\value{Object of class \code{"TotalVarIC"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{IC-class}}, \code{\link{ContIC}}}
\examples{
IC1 <- TotalVarIC()
plot(IC1)
}
\concept{influence curve}
\keyword{robust}
