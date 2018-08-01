\name{.checkEstClassForParamFamily-methods}
\docType{methods}
\alias{.checkEstClassForParamFamily}
\alias{.checkEstClassForParamFamily-methods}
\alias{.checkEstClassForParamFamily,GParetoFamily,Estimate-method}
\alias{.checkEstClassForParamFamily,GParetoFamily,LDEstimate-method}
\alias{.checkEstClassForParamFamily,GParetoFamily,MCEstimate-method}
\alias{.checkEstClassForParamFamily,GParetoFamily,kStepEstimate-method}
\alias{.checkEstClassForParamFamily,GParetoFamily,ORobEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamily,Estimate-method}
\alias{.checkEstClassForParamFamily,GEVFamily,MCEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamily,LDEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamily,kStepEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamily,ORobEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamilyMuUnknown,Estimate-method}
\alias{.checkEstClassForParamFamily,GEVFamilyMuUnknown,MCEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamilyMuUnknown,LDEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamilyMuUnknown,kStepEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamilyMuUnknown,ORobEstimate-method}
\title{ Methods for Function .checkEstClassForParamFamily in Package `RobExtremes' }
\description{.checkEstClassForParamFamily-methods}
%\usage{.checkEstClassForParamFamily(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GParetoFamily,Estimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GParetoFamily,MCEstimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GParetoFamily,LDEstimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GParetoFamily,kStepEstimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GEVFamily,Estimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GEVFamily,LDEstimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GEVFamily,MCEstimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GEVFamily,kStepEstimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GEVFamilyMuUnknown,Estimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GEVFamilyMuUnknown,LDEstimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GEVFamilyMuUnknown,MCEstimate}(PFam, estimator)
%\S4method{.checkEstClassForParamFamily}{GEVFamilyMuUnknown,kStepEstimate}(PFam, estimator)
%}
\arguments{
  \item{PFam}{a parametric family.}
  \item{estimator}{an estimator. }
  }
\details{
The respective methods can be used to cast an estimator to
a model-specific subclass with particular methods.
}

\value{
The \code{GParetoFamily,Estimate}-method returns the estimator cast to
S4 class \code{GPDEstimate},\cr
the \code{GParetoFamily,LDEstimate}-method cast to
S4 class  \code{GPDLDEstimate},\cr
the \code{GParetoFamily,MCEstimate}-method cast to
S4 class  \code{GPDMCEstimate},\cr
the \code{GParetoFamily,kStepEstimate}-method cast to
S4 class  \code{GPDkStepstimate},\cr
the \code{GEVFamily,Estimate}-method cast to
S4 class  \code{GEVEstimate},\cr
the \code{GEVFamily,LDEstimate}-method cast to
S4 class  \code{GEVLDEstimate},\cr
the \code{GEVFamily,MCEstimate}-method cast to
S4 class  \code{GEVMCEstimate},\cr
the \code{GEVFamily,kStepEstimate}-method cast to
S4 class  \code{GEVkStepstimate},\cr
the \code{GEVFamilyMuUnknown,Estimate}-method cast to
S4 class  \code{GEVEstimate},\cr
the \code{GEVFamilyMuUnknown,LDEstimate}-method cast to
S4 class  \code{GEVLDEstimate},\cr
the \code{GEVFamilyMuUnknown,MCEstimate}-method cast to
S4 class  \code{GEVMCEstimate},\cr
the \code{GEVFamilyMuUnknown,kStepEstimate}-method cast to
S4 class  \code{GEVkStepstimate}.\cr
}
\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
  }

\keyword{methods}

