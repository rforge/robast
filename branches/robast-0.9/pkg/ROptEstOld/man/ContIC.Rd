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
       lowerCase = NULL, neighborRadius = 0)
}
\arguments{
  \item{name}{ object of class \code{"character"}. }
  \item{CallL2Fam}{ object of class \code{"call"}:
      creates an object of the underlying L2-differentiable 
      parametric family. }
  \item{Curve}{ object of class \code{"EuclRandVarList"} }
  \item{Risks}{ object of class \code{"list"}: 
      list of risks; cf. \code{\link{RiskType-class}}. }
  \item{Infos}{ matrix of characters with two columns 
    named \code{method} and \code{message}: additional informations. }
  \item{clip}{ positive real: clipping bound. }
  \item{cent}{ real: centering constant }
  \item{stand}{ matrix: standardizing matrix }
  \item{lowerCase}{ optional constant for lower case solution. }
  \item{neighborRadius}{ radius of the corresponding (unconditional) 
    contamination neighborhood. }
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
\seealso{\code{\link{IC-class}}, \code{\link{ContIC}}}
\examples{
IC1 <- ContIC()
plot(IC1)
}
\concept{influence curve}
\keyword{robust}
