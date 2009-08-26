\name{ALEstimate-class}
\docType{class}
\alias{ALEstimate-class}
\alias{pIC}
\alias{pIC,ALEstimate-method}
\alias{asbias}
\alias{asbias,ALEstimate-method}
\alias{show,ALEstimate-method}
\alias{confint,ALEstimate,missing-method}
\alias{confint,ALEstimate,symmetricBias-method}
\alias{confint,ALEstimate,onesidedBias-method}
\alias{confint,ALEstimate,asymmetricBias-method}

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
Class \code{"Estimate"}, directly.
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
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distrMod]{Estimate-class}}}
\examples{
## prototype
new("ALEstimate")
}
\concept{estimate}
\keyword{classes}
