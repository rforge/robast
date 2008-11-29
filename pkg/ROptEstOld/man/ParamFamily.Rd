\name{ParamFamily}
\alias{ParamFamily}

\title{Generating function for ParamFamily-class}
\description{
  Generates an object of class \code{"ParamFamily"}.
}
\usage{
ParamFamily(name, distribution = Norm(), distrSymm, main = 0, 
            nuisance, trafo, param, props = character(0))
}
\arguments{
  \item{name}{ character string: name of family }
  \item{distribution}{ object of class \code{"Distribution"}: 
    member of the family }
  \item{distrSymm}{ object of class \code{"DistributionSymmetry"}: 
      symmetry of \code{distribution}. }
  \item{main}{ numeric vector: main parameter }
  \item{nuisance}{ numeric vector: nuisance parameter }
  \item{trafo}{ matrix: transformation of the parameters }
  \item{param}{ object of class \code{"ParamFamParameter"}:
    parameter of the family }
  \item{props}{ character vector: properties of the family }
}
\details{
  If \code{name} is missing, the default 
  \dQuote{"parametric family of probability measures"} is used. 
  In case \code{distrSymm} is missing it is set 
  to \code{NoSymmetry()}.
  If \code{param} is missing, the parameter is created via
  \code{main}, \code{nuisance} and \code{trafo} as described
  in \code{\link{ParamFamParameter}}.
}
\value{Object of class \code{"ParamFamily"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ParamFamily-class}}}
\examples{
F1 <- ParamFamily()
plot(F1)

## The function is currently defined as
function(name, distribution = Norm(), main = 0, nuisance, 
         trafo, param, props = character(0)){
    if(missing(name)) 
        name <- "parametric family of probability measures"
    if(missing(distrSymm)) distrSymm <- NoSymmetry()
    if(missing(param)) 
        param <- ParamFamParameter(name = paste("parameter of", name), 
                        main = main, nuisance = nuisance, trafo = trafo)
    return(new("ParamFamily", name = name, distribution = distribution, 
               distrSymm = distrSymm, param = param, props = props))
}
}
\concept{parametric family}
\keyword{distribution}
