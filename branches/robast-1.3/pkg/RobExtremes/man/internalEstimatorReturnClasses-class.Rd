\name{InternalEstimatorReturnClasses}
\docType{class}
\alias{InternalEstimatorReturnClasses-class}
\alias{GPDEstimate-class}
\alias{GEVEstimate-class}
\alias{GPDMCEstimate-class}
\alias{GEVMCEstimate-class}
\alias{GPDMDEstimate-class}
\alias{GEVMDEstimate-class}
\alias{GPDLDEstimate-class}
\alias{GEVLDEstimate-class}
\alias{GPDkStepEstimate-class}
\alias{GEVkStepEstimate-class}
\alias{GPDORobEstimate-class}
\alias{GEVORobEstimate-class}
\alias{GPDCvMMD.ALEstimate-class}
\alias{GEVCvMMD.ALEstimate-class}
\alias{GPDML.ALEstimate-class}
\alias{GEVML.ALEstimate-class}
\title{Internal Estimator Return Classes in 'RobExtremes'}
\description{S4 classes for return values of estimators
in package \pkg{RobExtremes} defined for internal
 purposes.}

 \section{Described classes}{
The S4 classes described here are \code{GPDEstimate}, \code{GEVEstimate},
 \code{GPDMCEstimate}, \code{GEVMCEstimate},
 \code{GPDMDEstimate}, \code{GEVMDEstimate},
 \code{GPDLDEstimate}, \code{GEVLDEstimate},
 \code{GPDkStepEstimate}, \code{GEVkStepEstimate},
 \code{GPDORobEstimate}, \code{GEVORobEstimate},
\code{GPDML.ALEstimate}, \code{GEVML.ALEstimate},
\code{GPDCvMMD.ALEstimate}, \code{GEVCvMMD.ALEstimate}.}

\section{Objects from the Class}{These classes are used internally
to provide specific S4 methods for different estimators later on;
thus, there are no generating functions.}
 \section{Slots}{
   All slots are inherited from parent classes.
}
\section{Extends}{
Classes \code{GPDEstimate}, \code{GEVEstimate} extend class \code{Estimate},
directly.\cr
Class \code{GPDMCEstimate} extends classes \code{GPDEstimate},
\code{MCEstimate}, directly.\cr
Class \code{GEVMCEstimate} extends classes \code{GEVEstimate},
\code{MCEstimate}, directly.\cr
Class \code{GPDMDEstimate} extends classes \code{GPDEstimate},
\code{MDEstimate}, directly.\cr
Class \code{GEVMDEstimate} extends classes \code{GEVEstimate},
\code{MDEstimate}, directly.\cr
Class \code{GPDMCALEstimate} extends classes \code{GPDEstimate},
\code{MCALEstimate}, directly.\cr
Class \code{GEVMCALEstimate} extends classes \code{GEVEstimate},
\code{MCALEstimate}, directly.\cr
Class \code{GPDLDEstimate} extends classes \code{GPDEstimate},
\code{LDEstimate}, directly.\cr
Class \code{GEVLDEstimate} extends classes \code{GEVEstimate},
\code{LDEstimate}, directly.\cr
Class \code{GPDkStepEstimate} extends classes \code{GPDEstimate},
\code{kStepEstimate}, directly.\cr
Class \code{GEVkStepEstimate} extends classes \code{GEVEstimate},
\code{kStepEstimate}, directly.\cr
Class \code{GPDORobEstimate} extends classes \code{GPDkStepEstimate},
\code{ORobEstimate}, directly.\cr
Class \code{GEVORobEstimate} extends classes \code{GEVkStepEstimate},
\code{ORobEstimate}, directly.\cr
Class \code{GPDML.ALEstimate} extends classes \code{GPDEstimate},
\code{ML.ALEstimate}, directly.\cr
Class \code{GEVML.ALEstimate} extends classes \code{GEVEstimate},
\code{ML.ALEstimate}, directly.\cr
Class \code{GPDCvMMD.ALEstimate} extends classes \code{GPDEstimate},
\code{CvMMD.ALEstimate}, directly.\cr
Class \code{GEVCvMMD.ALEstimate} extends classes \code{GEVEstimate},
\code{CvMMD.ALEstimate}, directly.\cr
}
%\references{}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link[distrMod]{Estimate-class}},
\code{\link[distrMod]{MCEstimate-class}},
\code{\link[RobAStBase]{kStepEstimate-class}},
\code{\link{LDEstimate-class}}}
\concept{GPD distribution}
\concept{GEV distribution}
\concept{S4 parameter class}
\keyword{models}
