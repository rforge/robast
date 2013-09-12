\name{movToRef-methods}
\docType{methods}
\alias{moveL2Fam2RefParam-methods}
\alias{moveICBackFromRefParam-methods}
\alias{moveICBackFromRefParam}
\alias{moveL2Fam2RefParam}
\alias{moveL2Fam2RefParam,L2ParamFamily-method}
\alias{moveL2Fam2RefParam,L2LocationFamily-method}
\alias{moveL2Fam2RefParam,L2ScaleFamily-method}
\alias{moveL2Fam2RefParam,L2LocationScaleFamily-method}
\alias{moveICBackFromRefParam,IC,L2ParamFamily-method}
\alias{moveICBackFromRefParam,IC,L2LocationFamily-method}
\alias{moveICBackFromRefParam,IC,L2ScaleFamily-method}
\alias{moveICBackFromRefParam,IC,L2LocationScaleFamily-method}
\alias{moveICBackFromRefParam,HampIC,L2ParamFamily-method}

\title{Methods for Functions moving from and to reference parameter in Package `ROptEst' }

\description{In \code{optIC} a gain in accuracy can be obtained when computing
 the optimally-robust ICs at a reference parameter of the model (instead of an
 arbtirary one). To this end, \code{moveL2Fam2RefParam} moved the model to
 the reference parameter and \code{moveICBackFromRefParam} moves the obtained
 optimal IC back to the original parameter.}

\usage{moveL2Fam2RefParam(L2Fam, ...)
       moveICBackFromRefParam(IC, L2Fam,...)
}

\arguments{
  \item{L2Fam}{object of class \code{L2ParamFamily}}
  \item{IC}{IC of class \code{HampIC}}
  \item{\dots}{further arguments to be passed on. }
}
\section{Methods}{\describe{
\item{moveL2Fam2RefParam}{\code{signature(L2Fam = "L2ParamFamily")}:
      returns \code{L2Fam} unchanged. }
\item{moveL2Fam2RefParam}{\code{signature(L2Fam = "L2LocationFamily")}:
      moves \code{L2Fam} to location \code{0}. }
\item{moveL2Fam2RefParam}{\code{signature(L2Fam = "L2ScaleFamily")}:
      moves \code{L2Fam} to location \code{0} and scale \code{1}. }
\item{moveL2Fam2RefParam}{\code{signature(L2Fam = "L2LocationScaleFamily")}:
      moves \code{L2Fam} to location \code{0} and scale \code{1}. }
\item{moveL2Fam2RefParam}{\code{signature(L2Fam = "L2LocationUnknownScaleFamily")}:
      moves \code{L2Fam} to location \code{0} and scale \code{1}. }
\item{moveL2Fam2RefParam}{\code{signature(L2Fam = "L2ScaleUnknownLocationFamily")}:
      moves \code{L2Fam} to location \code{0} and scale \code{1}. }
\item{moveICBackFromRefParam}{\code{signature(IC = "IC", L2Fam = "L2ParamFamily")}:
      returns \code{IC} unchanged. }
\item{moveICBackFromRefParam}{\code{signature(IC = "IC", L2Fam = "L2LocationFamily")}:
      moves IC in \code{IC} back to original location in \code{L2Fam}. }
\item{moveICBackFromRefParam}{\code{signature(IC = "IC", L2Fam = "L2ScaleFamily")}:
      moves IC in \code{IC} back to original location and scale in \code{L2Fam},
      rescaling risk where necessary. }
\item{moveICBackFromRefParam}{\code{signature(IC = "IC", L2Fam = "L2LocationScaleFamily")}:
      moves IC in \code{IC} back to original location and scale in \code{L2Fam},
      rescaling risk where necessary. }
\item{moveICBackFromRefParam}{\code{signature(IC = "IC", L2Fam = "L2LocationUnknownScaleFamily")}:
      moves IC in \code{IC} back to original location and scale in \code{L2Fam},
      rescaling risk where necessary. }
\item{moveICBackFromRefParam}{\code{signature(IC = "IC", L2Fam = "L2ScaleUnknownLocationFamily")}:
      moves IC in \code{IC} back to original location and scale in \code{L2Fam},
      rescaling risk where necessary. }
\item{moveICBackFromRefParam}{\code{signature(IC = "HampIC", L2Fam = "L2ParamFamily")}:
      moves IC in \code{IC} back to original location and scale in \code{L2Fam}
      (and in addition changes Lagrange multipliers accordingly),
      rescaling risk where necessary. }
}}
\value{
\item{\code{moveL2Fam2RefParam}}{the L2 Family transformed to reference parameter.}
\item{\code{moveICBackFromRefParam}}{the backtransformed IC.}
}

\details{\code{moveL2Fam2RefParam} and \code{moveICBackFromRefParam} are used
internally in functions \code{robest} and \code{roptest} to compute the
optimally robust influence function according to the arguments given to them.}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\seealso{\code{\link{robest}},\code{\link{optIC}}, \code{\link{radiusMinimaxIC}}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{classes}

