\name{InternalProbFamilyClasses-class}
\docType{class}
\alias{InternalProbFamilyClasses-class}
\alias{L2LocScaleShapeUnion-class}
\alias{ParamWithLocAndScaleAndShapeFamParameter-class}
\alias{locscaleshapename}
\alias{locscaleshapename-methods}
\alias{locscaleshapename,L2LocScaleShapeUnion-method}
\alias{locscalename,L2LocScaleShapeUnion-method}
\alias{scaleshapename,L2LocScaleShapeUnion-method}
\alias{locationname}
\alias{locationname-methods}
\alias{locationname,L2LocScaleShapeUnion-method}
\alias{scalename,L2LocScaleShapeUnion-method}
\alias{shapename}
\alias{shapename-methods}
\alias{shapename,L2LocScaleShapeUnion-method}
\alias{locscaleshapename<-}
\alias{locscaleshapename<--methods}
\alias{locscaleshapename<-,L2LocScaleShapeUnion-method}
\title{Internal Classes for Method Dispatch in 'ProbFamliy' Classes}
\description{Internal S4 classes for method dispatch in 'L2ParamFamily' and
'L2LocationFamily' (and friends) and in the respective parameter classes.}

\section{Described classes}{
In this file we describe classes \code{L2LocScaleShapeUnion} and
\code{ParamWithLocAndScaleAndShapeFamParameter}.\cr
Class \code{L2LocScaleShapeUnion} is a virtual class,
extending class \code{L2GroupParamFamily} with new slot \code{locscaleshapename}
(and, in fact, but not by S4 inheritance,
containing classes \code{L2ScaleShapeUnion} and \code{L2LocationScaleFamily}).
It is the parent class of class \code{GEVFamilyMuUnknown}.\cr
Class \code{ParamWithLocAndScaleAndShapeFamParameterUnion} is a virtual class
(union) containing classes \code{ParamWithScaleFamParameter} and
\code{ParamWithShapeFamParameter}.\cr
Class \code{ParamWithLocAndScaleAndShapeFamParameter} ``extends'' (no new
slots) class \code{ParamWithScaleAndShapeFamParameter}. It is the class
of the parameter in the class \code{GEVFamilyMuUnknown}.\cr
}
\section{Objects from these classes}{
  Objects are only generated internally by the mentioned generating functions. 
}
\section{Methods}{
  \describe{
    \item{locscaleshapename}{\code{signature(object = "L2LocationScaleShapeUnion")}:
      accesses the respective slot of the class }
    \item{locscalename}{\code{signature(object = "L2LocationScaleShapeUnion")}:
      accesses the location and scale part of the respective slot of the class }
    \item{scaleshapename}{\code{signature(object = "L2LocationScaleShapeUnion")}:
      accesses the scale and shape part of the respective slot of the class }
    \item{locationname}{\code{signature(object = "L2LocationScaleShapeUnion")}:
      accesses the location part of the respective slot of the class }
    \item{scalename}{\code{signature(object = "L2LocationScaleShapeUnion")}:
      accesses the scale part of the respective slot of the class }
    \item{shapename}{\code{signature(object = "L2LocationScaleShapeUnion")}:
      accesses the shape part of the respective slot of the class }
    \item{locscaleshapename<-}{\code{signature(object = "L2LocationScaleShapeUnion", value = "character")}:
      replaces the respective slot of the class }
}
}

\seealso{\code{\link[distrMod]{L2ParamFamily-class}},
\code{\link[distrMod:internalClassUnions-class]{L2ScaleShapeUnion-class}},
\code{\link[distrMod]{L2LocationScaleFamily-class}},
\code{\link[distrMod:ParamFamParameter-class]{ParamWithScaleAndShapeFamParameter-class}},
\code{\link[distrMod:ParamFamParameter-class]{ParamWithScaleFamParameter-class}},
\code{\link[distrMod:ParamFamParameter-class]{ParamWithShapeFamParameter-class}}.
}

\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\concept{parametric family}
\keyword{classes}
\keyword{models}
\keyword{internal}
