\name{updateNorm-methods}
\docType{methods}
\alias{updateNorm-methods}
\alias{updateNorm}
\alias{updateNorm,SelfNorm-method}
\title{ Methods for Function updateNorm in Package `ROptEst' }

\description{updateNorm-methods to update norm in IC-Algo}

\usage{updateNorm(normtype, ...)
\S4method{updateNorm}{SelfNorm}(normtype, L2, neighbor, biastype, Distr, V.comp, 
                                cent, stand,  w)     
}

\arguments{
  \item{normtype}{normtype of class \code{NormType}}
  \item{\dots}{ further arguments to be passed to specific methods.}
  \item{L2}{L2derivative}
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{biastype}{ object of class \code{"BiasType"} }
  \item{cent}{ optimal centering constant. }
  \item{stand}{ standardizing matrix. }
  \item{Distr}{ standardizing matrix. }
  \item{V.comp}{ matrix: indication which components of the standardizing
    matrix have to be computed. }
  \item{w}{object of class \code{RobWeight}; current weight}  
}
\section{Methods}{\describe{
\item{updateNorm}{\code{signature(normtype = "SelfNorm")}:
      udates the norm in the self-standardized case; just used
      internally in the opt-IC-Algorithm. }
}}
\value{
\item{updateNorm} an updated object of class \code{NormType}
}

\details{\code{updateNorm} is used internally in the opt-IC-algorithm to be
 able to work with a norm that depends on the current covariance 
 (\code{SelfNorm})}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
\seealso{\code{\link[distrMod]{NormType-class}}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{classes}

