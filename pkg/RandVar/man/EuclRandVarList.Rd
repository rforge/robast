\name{EuclRandVarList}
\alias{EuclRandVarList}

\title{Generating function for EuclRandVarList-class}
\description{
  Generates an object of class \code{"EuclRandVarList"}.
}
\usage{
EuclRandVarList(...)
}
\arguments{
  \item{\dots}{ Objects of class \code{"EuclRandVariable"} which
    shall form the list of Euclidean random variables. }
}
%\details{}
\value{Object of class \code{"EuclRandVarList"}}
% \references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{EuclRandVarList-class}}}
\examples{
L1 <- list(function(x){x}, function(x){x^2}, function(x){x^3}, function(x){x^4}, 
           function(x){x^5}, function(x){x^6})
L2 <- list(function(x){exp(x)}, function(x){abs(x)}, 
           function(x){sin(x)}, function(x){floor(x)})

R1 <- new("EuclRandVariable", Map = L2, Domain = Reals(), Range = Reals())
R2 <- EuclRandMatrix(Map = L1, ncol = 2, Domain = Reals(), dimension = 1)
R3 <- EuclRandMatrix(Map = L2, ncol = 2, Domain = Reals(), dimension = 1)

(RL1 <- EuclRandVarList(R1, R2, R3))
is(R1, "EuclRandVarList")
as(R1, "EuclRandVarList")
is(R2, "EuclRandVarList")
as(R2, "EuclRandVarList")

Map(exp(RL1)[[1]]) # "Math" group

## "Arith" group
Map((1 + RL1)[[1]])
Map((RL1 * 2)[[2]])
Map((RL1 / RL1)[[3]])

## The function is currently defined as
function(...){ 
    new("EuclRandVarList", list(...)) 
}
}
\concept{random variable}
\concept{random matrix}
\keyword{}
