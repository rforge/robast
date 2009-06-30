\name{RegTypeFamily}
\alias{RegTypeFamily}

\title{Generating function for RegTypeFamily-class}
\description{
  Generates an object of class \code{"RegTypeFamily"}.
}
\usage{
RegTypeFamily(name, distribution = LMCondDistribution(), distrSymm, 
        ErrorDistr = Norm(), ErrorSymm, main = 0, nuisance, trafo, 
        param, props = character(0), RegDistr = Norm(), RegSymm, 
        Regressor = RealRandVariable(c(function(x) {x}), Domain = Reals()))
}
\arguments{
  \item{name}{ name of the family }
  \item{distribution}{ conditional distribution (given the regressor) }
  \item{distrSymm}{ symmetry of \code{distribution} }
  \item{ErrorDistr}{ error distribution }
  \item{ErrorSymm}{ symmetry of \code{ErrorDistr} }
  \item{main}{ main parameter }
  \item{nuisance}{ optional nuisance parameter }
  \item{trafo}{ matrix: optional transformation of the parameter }
  \item{param}{ parameter of the family }
  \item{props}{ properties of the family }
  \item{RegDistr}{ regressor distribution }
  \item{RegSymm}{ symmetry of \code{RegDistr} }
  \item{Regressor}{ regressor }
}
\details{
  If \code{name} is missing, the default 
  \dQuote{regression type family} is used. 
  If \code{param} is missing, the parameter is created via
  \code{main}, \code{nuisance} and \code{trafo} as described
  in \code{ParamFamParameter}. In case \code{distrSymm}, 
  \code{ErrorSymm} or \code{RegSymm} is missing, they are set 
  to \code{NoSymmetry()}.
}
\value{Object of class \code{"RegTypeFamily"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{ParamFamily-class}}
\examples{
RegTypeFamily()
}
\keyword{models}
\keyword{robust}
