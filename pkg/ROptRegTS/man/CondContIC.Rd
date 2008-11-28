\name{CondContIC}
\alias{CondContIC}

\title{Generating function for CondContIC-class}
\description{
  Generates an object of class \code{"CondContIC"};
  i.e., an influence curves \eqn{\eta}{eta} of the form
  \deqn{\eta = (A\Lambda - a)\min(1,b/|A\Lambda - a|)}{eta = (A Lambda - a)min(1, b/|A Lambda - a|)}
  with clipping function \eqn{b}, centering function \eqn{a} and
  standardizing matrix \eqn{A}. \eqn{\Lambda}{Lambda} stands for
  the L2 derivative of the corresponding L2 differentiable 
  parametric family which can be created via \code{CallL2Fam}.
}
\usage{
CondContIC(name, CallL2Fam = call("L2RegTypeFamily"),
       Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x[1]*x[2]}), 
                                                Domain = EuclideanSpace(dimension = 2))), 
       Risks, Infos, 
       clip = RealRandVariable(Map = list(function(x){ Inf }), Domain = Reals()), 
       stand = as.matrix(1), 
       cent = EuclRandVarList(RealRandVariable(Map = list(function(x){numeric(length(x))}),
                                               Domain = EuclideanSpace(dimension = 2))), 
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
  \item{clip}{ object of class \code{"RealRandVariable"}: clipping function. }
  \item{cent}{ object of class \code{"EuclRandVarList"}: centering function. }
  \item{stand}{ matrix: standardizing matrix. }
  \item{lowerCase}{ optional constant for lower case solution. }
  \item{neighborRadius}{ radius of the corresponding conditional contamination 
    neighborhood. }
  \item{neighborRadiusCurve}{ radius curve of the corresponding conditional contamination 
    neighborhood. }
}
%\details{}
\value{Object of class \code{"CondContIC"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondIC-class}}, \code{\link{CondContIC-class}}}
\examples{
IC1 <- CondContIC()
IC1
}
\concept{influence curve}
\keyword{robust}

