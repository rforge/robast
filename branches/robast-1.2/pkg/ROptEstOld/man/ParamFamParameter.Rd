\name{ParamFamParameter}
\alias{ParamFamParameter}

\title{Generating function for ParamFamParameter-class}
\description{
  Generates an object of class \code{"ParamFamParameter"}.
}
\usage{
ParamFamParameter(name, main = numeric(0), nuisance, trafo)
}
\arguments{
  \item{name}{ character string: name of parameter }
  \item{main}{ numeric vector: main parameter }
  \item{nuisance}{ numeric vector: nuisance paramter }
  \item{trafo}{ matrix: transformation of the parameter }
}
\details{
  If \code{name} is missing, the default 
  \dQuote{"parameter of a parametric family of probability measures"}
  is used. If \code{nuisance} is missing, the nuisance parameter is 
  set to \code{NULL}. The number of columns of \code{trafo} have
  to be equal and the number of rows have to be not larger than 
  the sum of the lengths of \code{main} and \code{nuisance}. 
  If \code{trafo} is missing, no transformation to the parameter
  is applied; i.e., \code{trafo} is set to an identity matrix.
}
\value{Object of class \code{"ParamFamParameter"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ParamFamParameter-class}}}
\examples{
ParamFamParameter(main = 0, nuisance = 1, trafo = diag(c(1,2)))

## The function is currently defined as
function(name, main = numeric(0), nuisance, trafo){
    if(missing(name))
        name <- "parameter of a parametric family of probability measures"
    if(missing(nuisance))
        nuisance <- NULL
    if(missing(trafo))
        trafo <- diag(length(main)+length(nuisance))

    return(new("ParamFamParameter", name = name, main = main, 
               nuisance = nuisance, trafo = trafo))
}
}
\concept{parameter}
\keyword{robust}
