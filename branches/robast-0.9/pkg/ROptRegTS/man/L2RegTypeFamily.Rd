\name{L2RegTypeFamily}
\alias{L2RegTypeFamily}

\title{Generating function for L2RegTypeFamily-class}
\description{
  Generates an object of class \code{"RegTypeFamily"}.
}
\usage{
L2RegTypeFamily(name, distribution = LMCondDistribution(), distrSymm, 
        main = 0, nuisance, trafo, param, props = character(0), 
        L2deriv = EuclRandVarList(EuclRandVariable(
                  Map = list(function(x) {x[1] * x[2]}),
                  Domain = EuclideanSpace(dimension = 2),
                  dimension = 1)),
        ErrorDistr = Norm(), ErrorSymm, RegDistr = Norm(), RegSymm, 
        Regressor = RealRandVariable(Map = list(function(x) {x}), Domain = Reals()), 
        ErrorL2deriv = EuclRandVarList(RealRandVariable(Map = list(function(x) {x}), 
                                                        Domain = Reals())), 
        ErrorL2derivSymm, ErrorL2derivDistr, ErrorL2derivDistrSymm, FisherInfo)
}

\arguments{
  \item{name}{ name of the family }
  \item{distribution}{ conditional distribution (given the regressor) }
  \item{distrSymm}{ symmetry of \code{distribution} }
  \item{ErrorDistr}{ error distribution }
  \item{ErrorSymm}{ object of class \code{"DistributionSymmetry"}: 
    symmetry of \code{ErrorDistr} }
  \item{main}{ main parameter }
  \item{nuisance}{ optional nuisance parameter }
  \item{trafo}{ matrix: optional transformation of the parameter }
  \item{param}{ parameter of the family }
  \item{props}{ properties of the family }
  \item{RegDistr}{ regressor distribution }
  \item{RegSymm}{ object of class \code{"DistributionSymmetry"}: 
    symmetry of \code{RegDistr} }
  \item{Regressor}{ regressor }
  \item{L2deriv}{ object of class \code{"EuclRandVariable"}: L2 derivative }
  \item{ErrorL2deriv}{ object of class \code{"EuclRandVariable"}: 
    L2 derivative of \code{ErrorDistr} }
  \item{ErrorL2derivDistr}{ distribution of \code{ErrorL2deriv} }
  \item{ErrorL2derivSymm}{ object of class \code{"FunSymmList"}: 
    symmetry of \code{ErrorL2deriv} }
  \item{ErrorL2derivDistrSymm}{ object of class \code{"DistrSymmList"}: 
    symmetry of \code{ErrorL2derivDistr} }
  \item{FisherInfo}{ Fisher information matrix }
}
\details{
  If \code{name} is missing, the default
  \dQuote{L2 differentiable regression type family} is used.
  If \code{param} is missing, the parameter is created via
  \code{main}, \code{nuisance} and \code{trafo} as described
  in \code{ParamFamParameter}. In case \code{distrSymm}, 
  \code{ErrorSymm}, \code{RegSymm} is missing, they are
  set to \code{NoSymmetry()}. If \code{FisherInfo} is missing,
  it is computed via numerical integration.
}
\value{Object of class \code{"L2RegTypeFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2RegTypeFamily-class}}}
\examples{
L2RegTypeFamily()
}
\concept{regression}
\keyword{models}
\keyword{robust}
