\name{RobExtremes-package}
\alias{RobExtremes-package}
\alias{RobExtremes}
\docType{package}
\title{
RobExtremes -- optimally robust estimation for extreme value distributions
}
\description{
\pkg{RobExtremes} provides infrastructure for optimally robust estimation
for extreme value distributions, extending packages
\pkg{distr}, \pkg{distrEx}, \pkg{distrMod}, \pkg{robustbase},
\pkg{RobAStBase}, and \pkg{ROptEst}. Importing from packages \pkg{actuar},
\pkg{evd}, it provides S4 classes and methods for
\itemize{\item Gumbel distribution
         \item Generalized Extreme Value distribution (GEVD)
         \item Generalized Pareto distribution (GPD)
         \item Pareto distribution}
together with particular methods for expectations (as in package \pkg{distrEx})
---with particular integrators for GPD and GEVD---, variances, median,
IQR, skewness, kurtosis. In addition, extending estimators \code{Sn} and \code{Qn}
from package \pkg{robustbase}, we provide functionals for Sn and Qn. A new
asymmetric version of the \code{mad}, \code{kMAD} gives yet another robust
scale estimator (and functional). \cr

As to models, we provide the GPD model, together with (speeded-up) optimally
robust estimators, with LDEstimators (in general, and with \code{medkMAD},
\code{medSn} and \code{medQn} as particular ones) as starting estimators.
}

\details{
\tabular{ll}{
Package: \tab RobExtremes \cr
Version: \tab 0.8 \cr
Date: \tab 2012-05-17 \cr
Depends: \tab R (>= 2.14.0), methods, distr(>= 2.4), distrEx(>= 2.4), \cr
         \tab RandVar(>= 0.9), distrMod(>=2.5), RobAStBase(>=0.8-0),
         ROptEst(>=0.9), \cr
         \tab robustbase(>= 0.8-0), evd, actuar\cr
Imports: \tab RobAStRDA \cr
LazyLoad: \tab yes \cr
ByteCompile: \tab yes \cr
License: \tab LGPL-3 \cr
URL: \tab http://robast.r-forge.r-project.org/\cr
SVNRevision: \tab 450 \cr
}
}
\section{Classes}{
\preformatted{

[*]: there is a generating function with the same name in this package

##########################
Distribution Classes
##########################

"Distribution" (from distr)
|>"UnivariateDistribution" (from distr)
|>|>"AbscontDistribution" (from distr)
|>|>|>"Gumbel"    [*]
|>|>|>"Pareto"    [*]
|>|>|>"GPareto"   [*]
|>|>|>"GEVD"      [*]


##########################
Parameter Classes
##########################

"OptionalParameter" (from distr)
|>"Parameter" (from distr)
|>|>"GumbelParameter"
|>|>"ParetoParameter"
|>|>"GEVDParameter"
|>|>"GParetoParameter"

##########################
ProbFamily classes
##########################
slots: [<name>(<class>)]

"ProbFamily"                                  (from distrMod)
|>"ParamFamily"                               (from distrMod)
|>|>"L2ParamFamily"                           (from distrMod)
|>|>|>"L2GroupParamFamily"                    (from distrMod)
|>|>|>|>"ParetoFamily"                  [*]
|>|>|>|>"L2ScaleShapeUnion"                   (from distrMod)
|>|>|>|>|>"GParetoFamily"               [*]
|>|>|>|>|>"GEVFamily"                   [*]
|>|>|>|>|>"WeibullFamily"               [*]
|>|>|>|>"L2LocationScaleUnion"  /VIRTUAL/     (from distrMod)
|>|>|>|>|>"L2LocationFamily"                  (from distrMod)
|>|>|>|>|>|>"GumbelLocationFamily"      [*]
}
}

\section{Functions}{

\preformatted{
LDEstimator   Estimators for scale-shape models based on location and dispersion
medSn                    loc=median disp=Sn
medQn                    loc=median disp=Qn
medkMAD                  loc=median disp=kMAD
asvarMedkMAD               [asy. variance to MedkMADE]
PickandsEstimator        PickandsEstimator
asvarPickands              [asy. variance to PickandsE]
QuantileBCCEstimator     Quantile based estimator for the Weibull distribution
asvarQBCC                  [asy. variance to QuantileBCCE]
}}

\section{Generating Functions}{
\preformatted{

Distribution Classes
Gumbel                  Generating function for Gumbel-class
GEVD                    Generating function for GEVD-class
GPareto                 Generating function for GPareto-class
Pareto                  Generating function for Pareto-class

L2Param Families
ParetoFamily            Generating function for ParetoFamily-class
GParetoFamily           Generating function for GParetoFamily-class
GEVFamily               Generating function for GEVFamily-class
WeibullFamily           Generating function for WeibullFamily-class

}}
\section{Methods}{
\preformatted{

Functionals:
E                       Generic function for the computation of
                        (conditional) expectations
var                     Generic functions for the computation of
                        functionals
IQR                     Generic functions for the computation of
                        functionals
median                  Generic functions for the computation of
                        functionals
skewness                Generic functions for the computation of
                        functionals
kurtosis                Generic functions for the computation of
                        Functionals
Sn                       Generic function for the computation of
                        (conditional) expectations
Qn                     Generic functions for the computation of
                        functionals

}
}
\section{Constants}{
\preformatted{

EULERMASCHERONICONSTANT
APERYCONSTANT

}}

\section{Start-up-Banner}{
You may suppress the start-up banner/message completely by setting 
\code{options("StartupBanner"="off")} somewhere before loading this package by 
\code{library} or \code{require} in your R-code / R-session.

If option \code{"StartupBanner"} is not defined (default) or setting    
\code{options("StartupBanner"=NULL)} or 
\code{options("StartupBanner"="complete")} the complete start-up banner is 
displayed.

For any other value of option \code{"StartupBanner"} (i.e., not in 
\code{c(NULL,"off","complete")}) only the version information is displayed.

The same can be achieved by wrapping the \code{library} or \code{require}  call 
into either \code{suppressStartupMessages()} or 
\code{onlytypeStartupMessages(.,atypes="version")}.  
 
As for general \code{packageStartupMessage}'s, you may also suppress all
 the start-up banner by wrapping the \code{library} or \code{require} 
 call into \code{suppressPackageStartupMessages()} from 
 \pkg{startupmsg}-version 0.5 on.
 }

\section{Package versions}{
Note: The first two numbers of package versions do not necessarily reflect
 package-individual development, but rather are chosen for the
 RobAStXXX family as a whole in order to ease updating "depends"
 information.
}

\author{
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de},\cr
Matthias Kohl \email{Matthias.Kohl@stamats.de}, and \cr
Nataliya Horbenko \email{Nataliya.Horbenko@itwm.fraunhofer.de},\cr

\emph{Maintainer:}  Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
}
\references{
M. Kohl (2005): \emph{Numerical Contributions to the Asymptotic
Theory of Robustness.} PhD Thesis. Bayreuth. Available as
\url{http://www.stamats.de/ThesisMKohl.pdf}

P. Ruckdeschel, M. Kohl, T. Stabla, F. Camphausen (2006):
S4 Classes for Distributions, \emph{R News}, \emph{6}(2), 2-6. 
\url{http://CRAN.R-project.org/doc/Rnews/Rnews_2006-2.pdf}

M.~Kohl, P. Ruckdeschel, H.~Rieder (2010):
Infinitesimally Robust Estimation in General Smoothly Parametrized Models.
\emph{Stat. Methods Appl.}, \bold{19}, 333--354.\cr

Ruckdeschel, P. and Horbenko, N. (2011): Optimally-Robust Estimators in Generalized
Pareto Models. ArXiv 1005.1476. To appear at \emph{Statistics}.
DOI: 10.1080/02331888.2011.628022. \cr

Ruckdeschel, P. and Horbenko, N. (2012): Yet another breakdown point notion:
EFSBP --illustrated at scale-shape models. \emph{Metrika}, \bold{75}(8),
1025--1047.

%a more detailed manual for \pkg{distr}, \pkg{distrSim}, \pkg{distrTEst}, and \pkg{RobExtremes} may be downloaded from
%\url{http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf}\cr

a vignette for packages \pkg{distr}, \pkg{distrSim}, \pkg{distrTEst}, 
and \pkg{RobExtremes} is included into the mere documentation package \pkg{distrDoc}
and may be called by \code{require("distrDoc");vignette("distr")}

a homepage to this package is available under\cr
\url{http://robast.r-forge.r-project.org/}

}
\keyword{package}
\concept{S4 condition class}
\concept{S4 distribution class}
\concept{functional}
\concept{kurtosis}
\concept{median}
\concept{skewness}
\concept{kMAD}
\concept{IQR}
\concept{var}
\concept{E}
\concept{Sn}
\concept{Qn}
\concept{LDEstimator}
\concept{medkMAD}
\concept{medSn}
\concept{medQn}
\seealso{
\code{\link[distr:0distr-package]{distr-package}},
\code{\link[distrEx:0distrEx-package]{distrEx-package}},
\code{\link[distrMod:0distrMod-package]{distrMod-package}},
\code{\link[RobAStBase:0RobAStBase-package]{RobAStBase-package}},
\code{\link[ROptEst:0ROptEst-package]{ROptEst-package}}
}
