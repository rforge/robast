\name{InternalProbFamilyReturnClasses-class}
\docType{class}
\alias{InternalProbFamilyReturnClasses-class}
\alias{GumbelLocationFamily-class}
\alias{ParetoFamily-class}
\alias{GEVFamily-class}
\alias{GEVFamilyMuUnknown-class}
\alias{GParetoFamily-class}
\alias{WeibullFamily-class}

\title{Internal return classes for generating functions}
\description{Internal return classes for generating functions 'L2ParamFamily' and
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

\seealso{\code{\link[methods:BasicClasses]{numeric-class}},
\code{\link{L2ParamFamily-class}},
\code{\link{L2ScaleShapeUnion-class}},
}

\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\concept{parametric family}
\keyword{classes}
\keyword{models}
\keyword{internal}
