\name{.checkEstClassForParamFamily-methods}
\docType{methods}
\alias{.checkEstClassForParamFamily}
\alias{.checkEstClassForParamFamily-methods}
\alias{.checkEstClassForParamFamily,GParetoFamily,Estimate-method}
\alias{.checkEstClassForParamFamily,GParetoFamily,LDEstimate-method}
\alias{.checkEstClassForParamFamily,GParetoFamily,MCEstimate-method}
\alias{.checkEstClassForParamFamily,GParetoFamily,kStepEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamily,Estimate-method}
\alias{.checkEstClassForParamFamily,GEVFamily,MCEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamily,LDEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamily,kStepEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamilyMuUnknown,Estimate-method}
\alias{.checkEstClassForParamFamily,GEVFamilyMuUnknown,MCEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamilyMuUnknown,LDEstimate-method}
\alias{.checkEstClassForParamFamily,GEVFamilyMuUnknown,kStepEstimate-method}
\title{ Methods for Function .checkEstClassForParamFamily in Package `RobExtremes' }
\description{.checkEstClassForParamFamily-methods}
\usage{
.checkEstClassForParamFamily(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GParetoFamily,Estimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GParetoFamily,MCEstimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GParetoFamily,LDEstimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GParetoFamily,kStepEstimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GEVFamily,Estimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GEVFamily,LDEstimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GEVFamily,MCEstimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GEVFamily,kStepEstimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GEVFamilyMuUnknown,Estimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GEVFamilyMuUnknown,LDEstimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GEVFamilyMuUnknown,MCEstimate}(PFam, estimator)
\S4method{.checkEstClassForParamFamily}{GEVFamilyMuUnknown,kStepEstimate}(PFam, estimator)
}
\arguments{
  \item{PFam}{a parametric family.}
  \item{estimator}{an estimator. }
  }
\details{
The respective methods can be used to cast an estimator to
a model-specific subclass with particular methods.
}

\value{
The \code{GParetoFamily,Estimate}-method returns the estimator cast to \code{GPDEstimate},
the \code{GParetoFamily,LDEstimate}-method to \code{GPDLDEstimate},
the \code{GParetoFamily,MCEstimate}-method to \code{GPDMCEstimate},
the \code{GParetoFamily,kStepEstimate}-method to \code{GPDkStepstimate},
the \code{GEVFamily,Estimate}-method to \code{GEVEstimate},
the \code{GEVFamily,LDEstimate}-method to \code{GEVLDEstimate},
the \code{GEVFamily,MCEstimate}-method to \code{GEVMCEstimate},
the \code{GEVFamily,kStepEstimate}-method to \code{GEVkStepstimate}.
the \code{GEVFamilyMuUnknown,Estimate}-method to \code{GEVEstimate},
the \code{GEVFamilyMuUnknown,LDEstimate}-method to \code{GEVLDEstimate},
the \code{GEVFamilyMuUnknown,MCEstimate}-method to \code{GEVMCEstimate},
the \code{GEVFamilyMuUnknown,kStepEstimate}-method to \code{GEVkStepstimate}.
}
\author{
  Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
  }

\keyword{methods}

