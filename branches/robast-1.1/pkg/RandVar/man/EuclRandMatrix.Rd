\name{EuclRandMatrix}
\alias{EuclRandMatrix}

\title{Generating function for EuclRandMatrix-class}
\description{
  Generates an object of class \code{"EuclRandMatrix"}.
}
\usage{
EuclRandMatrix(Map = list(function(x){1}), nrow = 1, ncol = 1, 
            Domain = NULL, dimension = 1, Range)
}
\arguments{
  \item{Map}{ list of functions forming the map. }
  \item{nrow}{ number of rows.  }
  \item{ncol}{ number of columns. }
  \item{Domain}{ object of class \code{"OptionalrSpace"}: 
    domain of \code{Map} }
  \item{dimension}{ positive integer: dimension of the 
    range of \code{Map} }
  \item{Range}{ object of class \code{"OptionalrSpace"}: 
    range of \code{Map} }
}
%\details{}
\value{Object of class \code{"EuclRandMatrix"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{EuclRandMatrix-class}}}
\examples{
L1 <- list(function(x){x}, function(x){x^2}, function(x){x^3}, function(x){x^4}, 
           function(x){x^5}, function(x){x^6})
L2 <- list(function(x){exp(x)}, function(x){abs(x)}, 
           function(x){sin(x)}, function(x){floor(x)})

R1 <- EuclRandMatrix(Map = L1, nrow = 3, Domain = Reals(), dimension = 1)
R1[1:2, 2]
R1[1:2, 1:2]
Map(R1[1,2])
Map(t(R1)[2,1])

R2 <- EuclRandMatrix(Map = L2, ncol = 2, Domain = Reals(), dimension = 1)
(DL <- imageDistr(R2, Norm()))
plot(DL)

Map(gamma(R2)) # "Math" group

## "Arith" group
Map(2/R1)
Map(R2 * R2)


## The function is currently defined as
function(Map = list(function(x){1}), nrow = 1, ncol = 1,
                              Domain = NULL, dimension = 1) {
    if (missing(nrow)) 
        nrow <- ceiling(length(Map)/ncol)
    else if (missing(ncol)) 
        ncol <- ceiling(length(Map)/nrow)
    
    if(missing(Range))
        return(new("EuclRandMatrix", Map = Map, Domain = Domain, 
                   Range = EuclideanSpace(dimension = dimension),
                   Dim = as.integer(c(nrow, ncol))))
    else
        return(new("EuclRandMatrix", Map = Map, Domain = Domain, 
                   Range = Range, Dim = as.integer(c(nrow, ncol))))
}
}
\concept{random variable}
\concept{random matrix}
\keyword{classes}
