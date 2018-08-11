\name{ALEstimate-class}
\docType{class}
\alias{ALEstimate-class}
\alias{MCALEstimate-class}
\alias{ML.ALEstimate-class}
\alias{CvMMD.ALEstimate-class}
\alias{pIC}
\alias{pIC,ALEstimate-method}
\alias{pIC,MCEstimate-method}
\alias{pIC,MLEstimate-method}
\alias{pIC,CvMMDEstimate-method}
\alias{pIC,MCALEstimate-method}
\alias{pIC,ML.ALEstimate-method}
\alias{pIC,CvMMD.ALEstimate-method}
\alias{asbias}
\alias{asbias,ALEstimate-method}
\alias{show,ALEstimate-method}
\alias{show,MCALEstimate-method}
\alias{confint,ALEstimate,missing-method}
\alias{confint,ALEstimate,symmetricBias-method}
\alias{confint,ALEstimate,onesidedBias-method}
\alias{confint,ALEstimate,asymmetricBias-method}
\alias{.checkEstClassForParamFamily,ANY,MLEstimate-method}
\alias{.checkEstClassForParamFamily,ANY,CvMMDEstimate-method}
\alias{getPIC}
\alias{getPIC,ANY-method}
\alias{getPIC,MLEstimate-method}
\alias{getPIC,CvMMDEstimate-method}

\title{ALEstimate-class.}
\description{Class of asymptotically linear estimates.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ALEstimate", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}:
      name of the estimator. }
    \item{\code{estimate}}{Object of class \code{"ANY"}:
      estimate. }
    \item{\code{estimate.call}}{Object of class \code{"call"}:
      call by which estimate was produced.}
    \item{\code{samplesize}}{ object of class \code{"numeric"} ---
      the samplesize (only complete cases are counted)
      at which the estimate was evaluated. }
    \item{\code{completecases}}{ object of class \code{"logical"} ---
         complete cases at which the estimate was evaluated. }
    \item{\code{asvar}}{object of class \code{"OptionalNumericOrMatrix"}
      which may contain the asymptotic (co)variance of the estimator. }
    \item{\code{asbias}}{Optional object of class \code{"numeric"}:
      asymptotic bias. }
    \item{\code{pIC}}{Optional object of class \code{InfluenceCurve}:
      influence curve. }
    \item{\code{nuis.idx}}{ object of class \code{"OptionalNumeric"}: 
     indices of \code{estimate} belonging to the nuisance part. }
    \item{\code{fixed}}{ object of class \code{"OptionalNumeric"}:
     the fixed and known part of the parameter}
    \item{\code{Infos}}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
    \item{\code{trafo}}{ object of class \code{"list"}:
     a list with components \code{fct} and \code{mat} (see below). }
    \item{\code{untransformed.estimate}}{Object of class \code{"ANY"}:
      untransformed estimate.}
    \item{\code{untransformed.asvar}}{ object of class \code{"OptionalNumericOrMatrix"}
      which may contain the asymptotic (co)variance of the untransformed
      estimator. }
  }
}
\section{Extends}{
Class \code{ALEstimate} extends class \code{"Estimate"}, directly.
Class \code{MCALEstimate} extends classes
\code{"ALEstimate"}, and \code{"MCEstimate"} directly.
Class \code{ML.ALEstimate} extends classes
\code{"ALEstimate"}, and \code{"MLEstimate"} directly.
Class \code{CvM.ALEstimate} extends classes
\code{"ALEstimate"}, and \code{"CvMMDEstimate"} directly.
The last two classes are to be used for method dispatch, later;
they have an identical slot structure to class \code{MCALEstimate}.
}
\section{Methods}{
  \describe{
    \item{pIC}{\code{signature(object = "ALEstimate")}: 
      accessor function for slot \code{pIC}. }

    \item{show}{\code{signature(object = "ALEstimate")} }

    \item{confint}{\code{signature(object = "ALEstimate", method = "missing")}: 
      compute asymptotic (LAN-based) confidence interval neglecting any bias. }

    \item{confint}{\code{signature(object = "ALEstimate", method = "symmetricBias")}: 
      compute asymptotic (LAN-based) confidence interval incorporating bias
      symmetrically. }

    \item{confint}{\code{signature(object = "ALEstimate", method = "onesidedBias")}: 
      compute asymptotic (LAN-based) confidence interval incorporating bias
      one-sided; i.e., positive or negative, respectively. }

    \item{confint}{\code{signature(object = "ALEstimate", method = "asymmetricBias")}: 
      compute asymptotic (LAN-based) confidence interval incorporating bias
      asymmetrically. }
  }
}
\details{The (return value) class of an estimator is of class \code{ALEstimate}
   if it is asymptotically linear; then it has an influence function
   (implemented in slot \code{pIC}) and so all the diagnostics for influence
   functions are available; in addition it is asymptotically normal, so
   we can (easily) deduce asymptotic covariances, hence may use these
   in confidence intervals; in particular, the return values of \code{kStepEstimator}
   \code{oneStepEstimator} (and \code{roptest}, \code{robest}, \code{RMXEstimator},
   \code{MBREstimator}, \code{OBREstimator}, \code{OMSEstimator} in package
   'ROptEst') are objects of (subclasses of) this class.

   As the return value of \code{CvMMDEEstimator} (or \code{MDEstimator} with
   \code{CvMDist} or \code{CvMDist2} as distance) is asymptotically linear,
   there is class \code{MCALEstimate} extending \code{MCEstimate} by
   extra slots \code{pIC} and \code{asbias} (only filled optionally with
   non-\code{NULL} values). Again all the diagnostics for influence
   functions are then available. Classes \code{ML.ALEstimate} and
   class \code{CvMMD.ALEstimate} are nominal subclasses of class
   \code{MCALEstimate}, nominal in the sense that they have no extra slots,
   but they might have particular methods later on.

   Helper method \code{getPIC} by means of the estimator class, and, in
   case of estimators of class \code{CvMMDEstimate}, also the name
   (in slot \code{name}) produces the (partial) influence function:
    calling \code{.CvMMDCovariance} -- either directly or through wrapper
   \code{.CvMMDCovarianceWithMux}. This is used  in the corresponding
   \code{.checkEstClassForParamFamily} method, which coerces \code{object} from
   class \code{"MCEstimate"} to \code{"MCALEstimate"}.
   }

%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de} and
Peter Ruckdeschel \email{Peter.Ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link[distrMod]{Estimate-class}}}
\examples{
## prototype
new("ALEstimate")

## data example
set.seed(123)
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

mle <- MLEstimator(x,G)
(picM <- pIC(mle))

## Kolmogorov(-Smirnov) minimum distance estimator
ke <- KolmogorovMDEstimator(x = x, ParamFamily = G)
pIC(ke) ## gives NULL

## von Mises minimum distance estimator with default mu

\donttest{ ## to save time for CRAN
system.time(me <- CvMMDEstimator(x = x, ParamFamily = G))
str(me@pIC) ## a call
system.time(pIC0 <- pIC(me))
str(me@pIC) ## now filled
}


}
\concept{estimate}
\keyword{classes}
