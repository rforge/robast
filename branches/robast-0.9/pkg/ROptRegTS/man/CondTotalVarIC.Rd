\name{CondTotalVarIC}
\alias{CondTotalVarIC}

\title{Generating function for CondTotalVarIC-class}
\description{
  Generates an object of class \code{"CondTotalVarIC"};
  i.e., an influence curves \eqn{\eta}{eta} of the form
  \deqn{\eta = max(c(x), min(A x\Lambda_f, b(x)))}{
  eta = max(c(x), min(A x Lambda_f, b(x)))}
  with lower clipping function \eqn{c}, upper clipping function \eqn{b} and
  standardizing matrix \eqn{A}. \eqn{\Lambda_f}{Lambda_f} stands for
  the L2 derivative of the corresponding error distribution.
}
\usage{
CondTotalVarIC(name, CallL2Fam = call("L2RegTypeFamily"), 
        Curve = EuclRandVarList(RealRandVariable(
                Map = list(function(x) {x[1] * x[2]}),
                Domain = EuclideanSpace(dimension = 2))),
        Risks, Infos, 
        clipUp = RealRandVariable(Map = list(function(x) {Inf}), Domain = Reals()), 
        stand = as.matrix(1), 
        clipLo = RealRandVariable(Map = list(function(x) {-Inf}), Domain = Reals()), 
        lowerCase = NULL, neighborRadius = 0, neighborRadiusCurve = function(x){1})
}
\arguments{
  \item{name}{ object of class \code{"character"}. }
  \item{CallL2Fam}{ object of class \code{"call"}:
      creates an object of the underlying L2-differentiable 
      regression type family. }
  \item{Curve}{ object of class \code{"EuclRandVarList"} }
  \item{Risks}{ object of class \code{"list"}: 
      list of risks; cf. \code{RiskType-class}. }
  \item{Infos}{ matrix of characters with two columns 
    named \code{method} and \code{message}: additional informations. }
  \item{clipUp}{ object of class \code{"RealRandVariable"}:
    upper clipping function. }
  \item{clipLo}{ object of class \code{"RealRandVariable"}:
    lower clipping function. }
  \item{stand}{ matrix: standardizing matrix. }
  \item{lowerCase}{ optional constant for lower case solution. }
  \item{neighborRadius}{ radius of the corresponding conditional
    total variation neighborhood. }
  \item{neighborRadiusCurve}{ radius curve of the corresponding conditional
    total variation neighborhood. }
}
%\details{}
\value{Object of class \code{"CondTotalVarIC"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondIC-class}}, \code{\link{CondTotalVarIC-class}}}
\examples{
IC1 <- CondTotalVarIC()
IC1
}
\concept{influence curve}
\keyword{robust}

