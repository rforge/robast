\name{Av1CondTotalVarIC}
\alias{Av1CondTotalVarIC}

\title{Generating function for Av1CondTotalVarIC-class}
\description{
  Generates an object of class \code{"Av1CondContIC"};
  i.e., an influence curves \eqn{\eta}{eta} of the form
  \deqn{\eta = A x \Lambda_f\min(1,\max(c(x)/(|Ax|\Lambda_f), (c(x)+b)/(|Ax|\Lambda_f)))}{
  eta = A x Lambda_f min(1, max(c(x)/(|Ax|Lambda_f), (c(x) + b)/(|Ax|Lambda_f)))}
  with lower clipping function \eqn{c}, standardized bias \eqn{b} and
  standardizing matrix \eqn{A}. \eqn{\Lambda_f}{Lambda_f} stands for
  the L2 derivative of the corresponding error distribution.
}
\usage{
Av1CondTotalVarIC(name, CallL2Fam = call("L2RegTypeFamily"), 
        Curve = EuclRandVarList(RealRandVariable(Map = list(function(x) {x[1] * x[2]}), 
                                                 Domain = EuclideanSpace(dimension = 2))), 
        Risks, Infos, clipUp = Inf, stand = as.matrix(1), 
        clipLo = RealRandVariable(Map = list(function(x) {-Inf}), 
                                  Domain = EuclideanSpace(dimension = 1)), 
                                  lowerCase = NULL, neighborRadius = 0)
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
  \item{clipUp}{ positive real: standardized bias. }
  \item{clipLo}{ object of class \code{"RealRandVariable"}:
    lower clipping function. }
  \item{stand}{ matrix: standardizing matrix. }
  \item{lowerCase}{ optional constant for lower case solution. }
  \item{neighborRadius}{ radius of the corresponding (unconditional) 
    contamination neighborhood. }
}
%\details{}
\value{Object of class \code{"Av1CondTotalVarIC"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondIC-class}}, \code{\link{Av1CondTotalVarIC-class}}}
\examples{
IC1 <- Av1CondTotalVarIC()
IC1
}
\concept{influence curve}
\keyword{robust}

