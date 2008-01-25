\name{Av2CondContIC}
\alias{Av2CondContIC}

\title{Generating function for Av2CondContIC-class}
\description{
  Generates an object of class \code{"Av2CondContIC"};
  i.e., an influence curves \eqn{\eta}{eta} of the form
  \deqn{\eta = A K^{-1} x(\Lambda_f - z)\min(1,c/|\Lambda_f - z|)}{eta = A K^(-1) x(Lambda_f - z)min(1, c/|Lambda_f - a|)}
  with \eqn{K=E xx^\tau}{K=E x t(x)}, clipping bound \eqn{c}, centering constant 
  \eqn{z} and standardizing constant \eqn{A}. \eqn{\Lambda_f}{Lambda_f} stands for
  the L2 derivative of the corresponding error distribution.
}
\usage{
Av2CondContIC(name, CallL2Fam = call("L2RegTypeFamily"), 
        Curve = EuclRandVarList(RealRandVariable(Map = list(function(x) {x[1] * x[2]}),    
                                                 Domain = EuclideanSpace(dimension = 2))), 
        Risks, Infos, clip = Inf, stand = 1, cent = 0, lowerCase = NULL, neighborRadius = 0)
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
  \item{clip}{ positive real: clipping bound. }
  \item{cent}{ real: centering constant }
  \item{stand}{ real: standardizing constant }
  \item{lowerCase}{ optional constant for lower case solution. }
  \item{neighborRadius}{ radius of the corresponding (unconditional) 
    contamination neighborhood. }
}
%\details{}
\value{Object of class \code{"Av2CondContIC"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondIC-class}}, \code{\link{Av2CondContIC-class}}}
\examples{
IC1 <- Av2CondContIC()
IC1
}
\concept{influence curve}
\keyword{}
