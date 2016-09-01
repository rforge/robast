\name{L2ParamFamily}
\alias{L2ParamFamily}

\title{Generating function for L2ParamFamily-class}
\description{
  Generates an object of class \code{"L2ParamFamily"}.
}
\usage{
L2ParamFamily(name, distribution = Norm(), distrSymm,
        main = 0, nuisance, trafo, param, props = character(0), 
        L2deriv = EuclRandVarList(RealRandVariable(list(function(x) {x}), 
                                                   Domain = Reals())), 
        L2derivSymm, L2derivDistr, L2derivDistrSymm, FisherInfo)
}
\arguments{
  \item{name}{ character string: name of the family }
  \item{distribution}{ object of class \code{"Distribution"}:
    member of the family }
  \item{distrSymm}{ object of class \code{"DistributionSymmetry"}: 
      symmetry of \code{distribution}. }
  \item{main}{ numeric vector: main parameter }
  \item{nuisance}{ numeric vector: nuisance parameter }
  \item{trafo}{ matrix: transformation of the parameter }
  \item{param}{ object of class \code{"ParamFamParameter"}:
    parameter of the family }
  \item{props}{ character vector: properties of the family }
  \item{L2deriv}{ object of class \code{"EuclRandVariable"}: 
    L2 derivative of the family }
  \item{L2derivSymm}{ object of class \code{"FunSymmList"}: 
    symmetry of the maps contained in \code{L2deriv} }
  \item{L2derivDistr}{ object of class \code{"UnivarDistrList"}: 
    distribution of \code{L2deriv} }
  \item{L2derivDistrSymm}{ object of class \code{"DistrSymmList"}: 
    symmetry of the distributions contained in \code{L2derivDistr} }
  \item{FisherInfo}{ object of class \code{"PosDefSymmMatrix"}:
    Fisher information of the family }
}
\details{
  If \code{name} is missing, the default 
  \dQuote{L2 differentiable parametric family of probability measures} 
  is used. In case \code{distrSymm} is missing it is set to
  \code{NoSymmetry()}.
  If \code{param} is missing, the parameter is created via
  \code{main}, \code{nuisance} and \code{trafo} as described
  in \code{\link{ParamFamParameter}}. In case \code{L2derivSymm} is
  missing, it is filled with an object of class \code{FunSymmList} 
  with entries \code{NonSymmetric()}. In case \code{L2derivDistr} is missing, 
  it is computed via \code{imageDistr}. If \code{L2derivDistrSymm} is missing, 
  it is set to an object of class \code{DistrSymmList} with entries 
  \code{NoSymmetry()}. In case \code{FisherInfo} is missing, it is computed 
  from \code{L2deriv} using \code{E}.
}
\value{Object of class \code{"L2ParamFamily"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}}
\examples{
F1 <- L2ParamFamily()
plot(F1)
}
\concept{parametric family}
\keyword{models}
