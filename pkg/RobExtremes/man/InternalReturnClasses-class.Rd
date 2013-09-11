\name{InternalReturnClasses-class}
\docType{class}
\alias{InternalReturnClasses-class}
\alias{GumbelLocationFamily-class}
\alias{ParetoFamily-class}
\alias{GEVFamily-class}
\alias{GEVFamilyMuUnknown-class}
\alias{GParetoFamily-class}
\alias{WeibullFamily-class}

\title{Internal return classes for generating functions}
\description{internal return classes for generating functions 'L2ParamFamily' and 
'L2LocationFamily' (and friends); used for particular method dispatch only}

\section{Described classes}{
In this file we describe classes \code{GParetoFamily}, \code{GEVFamily},
\code{GEVFamilyMuUnknown}, \code{WeibullFamily} all ``extending''
(no new slots!) class union
\code{"L2LocationScaleShapeUnion"} and \code{ParetoFamily} ``extending''
 (no new slots!) class \code{L2ParamFamily}.
}
\section{Objects from these classes}{
  Objects are only generated internally by the mentioned generating functions. 
}

\section{Methods}{
   not yet done...
}
\seealso{\code{\link[methods:BasicClasses]{numeric-class}},
\code{\link{L2ParamFamily-class}},
\code{\link{L2ScaleShapeUnion-class}},
}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\concept{parametric family}
\keyword{classes}
\keyword{models}
\keyword{internal}
