\name{movToRef-methods}
\docType{methods}
\alias{moveL2Fam2RefParam}
\alias{moveICBackFromRefParam}
\alias{moveL2Fam2RefParam-methods}
\alias{moveICBackFromRefParam-methods}
\alias{moveL2Fam2RefParam,L2ScaleShapeUnion-method}
\alias{moveICBackFromRefParam,IC,L2ScaleShapeUnion-method}
\alias{moveICBackFromRefParam,IC,L2LocScaleShapeUnion-method}

\title{Methods for Functions moving from and to reference parameter in Package
    `RobExtremes' }

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
  \item{\dots}{further arguments to be passed to particular methods}
}
\section{Methods}{\describe{
\item{moveL2Fam2RefParam}{\code{signature(L2Fam = "L2ScaleShapeUnion")}:
      moves \code{L2Fam} to scale \code{1} (and, if existing location to
       \code{0}). }
\item{moveICBackFromRefParam}{\code{signature(IC = "IC", L2Fam = "L2ScaleShapeUnion")}:
      moves IC in \code{IC} back to original location and scale in \code{L2Fam}
      (and in addition changes Lagrange multipliers accordingly),
      rescaling risk where necessary. }
\item{moveICBackFromRefParam}{\code{signature(IC = "IC", L2Fam = "L2LocScaleShapeUnion")}:
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
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
\seealso{\code{\link{robest}},\code{\link{optIC}}, \code{\link{radiusMinimaxIC}}}
%\examples{}
\concept{asymptotic risk}
\concept{risk}
\keyword{classes}

