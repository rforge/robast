\name{RealRandVariable}
\alias{RealRandVariable}

\title{Generating function for RealRandVariable-class}
\description{
  Generates an object of class \code{"RealRandVariable"}.
}
\usage{
RealRandVariable(Map = list(function(x) {1}), Domain = NULL, Range)
}
\arguments{
  \item{Map}{ list of functions forming the map. }
  \item{Domain}{ domain of \code{Map}: object of 
    class \code{"OptionalrSpace"}. }
  \item{Range}{ range of \code{Map}: object of 
    class \code{"Reals"}. }
}
%\details{}
\value{Object of class \code{"RealRandVariable"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{RealRandVariable-class}}}
\examples{
RealRandVariable(Map = list(function(x){x}), Domain = Reals())

## The function is currently defined as
function(Map = list(function(x){1}), Domain = NULL, Range) {
    if(missing(Range)) Range <- Reals()
    if(!is(Range, "Reals"))
        stop("'Range' has to be of class 'Reals'")

    return(new("RealRandVariable", Map = Map, 
               Domain = Domain, Range = Reals()))
}
}
\concept{random variable}
\keyword{classes}
