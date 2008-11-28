\name{generateIC-methods}
\docType{methods}
\alias{generateIC-methods}
\alias{generateIC,ContNeighborhood,L2RegTypeFamily-method}
\alias{generateIC,TotalVarNeighborhood,L2RegTypeFamily-method}

\title{Methods for Function generateIC in Package `ROptRegTS'}
\description{
 Methods for function \code{generateIC} in package \pkg{ROptRegTS}.
}
\section{Methods}{
\describe{
  \item{neighbor = "ContNeighborhood", L2Fam = "L2RegTypeFamily"}{ 
    generate an object of class \code{"ContIC"}. Rarely called directly. }

  \item{neighbor = "TotalVarNeighborhood", L2Fam = "L2RegTypeFamily"}{ 
    generate an object of class \code{"TotalVarIC"}. Rarely called directly. }
}}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\seealso{\code{\link[ROptEst]{generateIC}}}
\concept{influence curve}
\keyword{methods}
\keyword{robust}
