\name{ORobEstimate-class}
\docType{class}
\alias{ORobEstimate-class}
\alias{steps}
\alias{steps,ORobEstimate-method}
\alias{ksteps}
\alias{ksteps,ORobEstimate-method}
\alias{uksteps}
\alias{uksteps,ORobEstimate-method}
\alias{ICList}
\alias{ICList,ORobEstimate-method}
\alias{pICList}
\alias{pICList,ORobEstimate-method}
\alias{roptestCall}
\alias{roptestCall,ORobEstimate-method}
\alias{robestCall}
\alias{robestCall,ORobEstimate-method}
\alias{start,ORobEstimate-method}
\alias{startval}
\alias{startval,ORobEstimate-method}
\alias{ustartval}
\alias{ustartval,ORobEstimate-method}
\alias{show,ORobEstimate-method}

\title{ORobEstimate-class.}
\description{Class of optimally robust asymptotically linear estimates.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ORobEstimate", ...)}.
  More frequently they are created as results of functions
  \code{roptest}, \code{MBREstimator}, \code{RMXEstimator}, or
  \code{OMSEstimator}.
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}:
      name of the estimator. [*]}
    \item{\code{estimate}}{Object of class \code{"ANY"}:
      estimate. [*]}
    \item{\code{estimate.call}}{Object of class \code{"call"}:
      call by which estimate was produced. [*]}
    \item{\code{samplesize}}{ object of class \code{"numeric"} ---
      the samplesize (only complete cases are counted)
      at which the estimate was evaluated. [*]}
    \item{\code{completecases}:}{ object of class \code{"logical"} ---
         complete cases at which the estimate was evaluated. [*]}
    \item{\code{asvar}}{ object of class \code{"OptionalNumericOrMatrix"}
      which may contain the asymptotic (co)variance of the estimator. [*]}
    \item{\code{asbias}}{Optional object of class \code{"numeric"}:
      asymptotic bias. [*]}
    \item{\code{pIC}}{Optional object of class \code{InfluenceCurve}:
      influence curve. [*]}
    \item{\code{nuis.idx}}{ object of class \code{"OptionalNumeric"}: 
     indices of \code{estimate} belonging to the nuisance part. [*]}
    \item{\code{fixed}}{ object of class \code{"OptionalNumeric"}:
     the fixed and known part of the parameter. [*]}
    \item{\code{steps}}{Object of class \code{"integer"}: number
      of steps. [*]}
    \item{\code{Infos}}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. [*]}
    \item{\code{trafo}}{ object of class \code{"list"}:
     a list with components \code{fct} and \code{mat} (see below). [*]}
    \item{\code{untransformed.estimate}:}{Object of class \code{"ANY"}:
      untransformed estimate. [*]}
    \item{\code{untransformed.asvar}:}{ object of class \code{"OptionalNumericOrMatrix"}
      which may contain the asymptotic (co)variance of the untransformed
      estimator. [*]}
    \item{\code{pICList}}{Optional object of class \code{"OptionalpICList"}:
      the list of (intermediate) (partial) influence curves used;
      only filled when called from \code{ORobEstimator} with argument
      \code{withPICList==TRUE}. [*]}
    \item{\code{ICList}}{Optional object of class \code{"OptionalpICList"}:
      the list of (intermediate) (total) influence curves used;
      only filled when called from \code{ORobEstimator} with argument
      \code{withICList==TRUE}. [*]}
    \item{\code{start}}{The argument \code{start} --- of class \code{"StartClass"}
    used in call to \code{ORobEstimator}. [*]}
    \item{\code{startval}}{Object of class \code{matrix}:
      the starting value with which the k-step Estimator was initialized
      (in \eqn{p}-space / transformed). [*]}
    \item{\code{ustartval}}{Object of class \code{matrix}:
      the starting value with which the k-step Estimator was initialized
      (in \eqn{k}-space / untransformed). [*]}
    \item{\code{ksteps}}{Object of class \code{"OptionalMatrix"}:
    the intermediate estimates (in \eqn{p}-space) for the parameter;
    only filled when called from \code{ORobEstimator}. [*]}
    \item{\code{uksteps}}{Object of class \code{"OptionalMatrix"}:
    the intermediate estimates (in \eqn{k}-space) for the parameter;
    only filled when called from \code{ORobEstimator}. [*]}
    \item{\code{robestcall}}{Object of class \code{"OptionalCall"}, i.e.,
     a \code{call} or \code{NULL}: only filled when called from \code{roptest}. [*]}
    \item{\code{roptestcall}}{Object of class \code{"OptionalCall"}, i.e.,
     a \code{call} or \code{NULL}: only filled when called from \code{roptest},
     \code{MBREstimator}, \code{RMXEstimator}, or \code{OMSEstimator}.}
  }
}
\section{Extends}{
Class \code{"kStepEstimate"}, directly.\cr
Class \code{"ALEstimate"} and class \code{"Estimate"}, by
class \code{"kStepstimate"}. All slots and methods marked with [*] are inherited.
}
\section{Methods}{
  \describe{
    \item{steps}{\code{signature(object = "ORobEstimate")}:
      accessor function for slot \code{steps}. [*]}

    \item{ksteps}{\code{signature(object = "ORobEstimate")}:
      accessor function for slot \code{ksteps}; has additional argument
      \code{diff}, defaulting to \code{FALSE}; if the latter is \code{TRUE},
      the starting value from slot \code{startval} is prepended as first column;
      otherwise we return the corresponding increments in each step. [*]}

    \item{uksteps}{\code{signature(object = "ORobEstimate")}:
      accessor function for slot \code{uksteps}; has additional argument
      \code{diff}, defaulting to \code{FALSE}; if the latter is \code{TRUE},
      the starting value from slot \code{ustartval} is prepended as first column;
      otherwise we return the corresponding increments in each step. [*]}

    \item{start}{\code{signature(object = "ORobEstimate")}:
      accessor function for slot \code{start}. [*]}

    \item{startval}{\code{signature(object = "ORobEstimate")}:
      accessor function for slot \code{startval}. [*]}

    \item{ustartval}{\code{signature(object = "ORobEstimate")}:
      accessor function for slot \code{startval}. [*]}

    \item{ICList}{\code{signature(object = "ORobEstimate")}:
      accessor function for slot \code{ICList}. [*]}

    \item{pICList}{\code{signature(object = "ORobEstimate")}:
      accessor function for slot \code{pICList}. [*]}

    \item{robestCall}{\code{signature(object = "ORobEstimate")}:
      accessor function for slot \code{robestCall}. [*]}

    \item{roptestCall}{\code{signature(object = "ORobEstimate")}:
      accessor function for slot \code{roptestCall}. }

    \item{show}{\code{signature(object = "ORobEstimate")}: a show method; [*]}
  }
}
%\references{}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{ALEstimate-class}}, \code{\link{kStepEstimate-class}}}
%\examples{}
\concept{estimate}
\keyword{classes}
