\name{RobExtremes-package}
\alias{RobExtremes-package}
\alias{RobExtremes}
\docType{package}
\title{
RobExtremes -- Optimally Robust Estimation for Extreme Value Distributions
}
\description{
\pkg{RobExtremes} provides infrastructure for speeded-up optimally robust
estimation (i.e., MBRE, OMSE, RMXE) for extreme value distributions,
extending packages \pkg{distr}, \pkg{distrEx}, \pkg{distrMod},
\pkg{robustbase}, \pkg{RobAStBase}, and \pkg{ROptEst}.}

\section{Distributions}{

Importing from packages \pkg{actuar},
\pkg{evd}, it provides S4 classes and methods for the
\itemize{\item Gumbel distribution
         \item Generalized Extreme Value distribution (GEVD)
         \item Generalized Pareto distribution (GPD)
         \item Pareto distribution}
%\cr
}
\section{Functionals for Distributions}{

These distributions come together with particular methods for expectations.
I.e., a functional E() as in package \pkg{distrEx}, which as first argument
takes the distribution, and, optionally, can take as second argument
a function which then is used as integrand. These particular methods are
available for the GPD, Pareto, Gamma, Weibull, and GEV disdribution and use
integration on the quantile scale, i.e.,
\deqn{\mathop{E}[X]=\int_0^1 q^X(s)\,ds}{E[X] = integral_0^1 q^X(s) ds}
where \eqn{q^X}{q^X} is the quantile function of X.
In addition, where they exist, we provide closed from expressions for
variances, median, IQR, skewness, kurtosis. \cr

In addition, extending estimators \code{Sn} and \code{Qn} from package
\pkg{robustbase}, we provide functionals for Sn and Qn. A new
asymmetric version of the \code{mad}, \code{kMAD} gives yet another robust
scale estimator (and functional). %\cr
}

\section{Models and Estimators}{

As to models, we provide the
\itemize{
\item GPD model (with known threshold), together with (speeded-up) optimally
robust estimators, with LDEstimators (in general, and with \code{medkMAD},
\code{medSn} and \code{medQn} as particular ones) and Pickands' estimator as
starting estimators.
\item GEVD model (with known or unknown threshold), together with (speeded-up)
optimally robust estimators, with LDEstimators (see above) and Pickands'
estimator as starting estimators.
\item Pareto model
\item Weibull model
\item Gamma model
}
and for each of these, we provide  speeded-up optimally robust estimation
(i.e., MBRE, OMSE, RMXE).\cr

We robust (high-breakdown) starting estimators for
\itemize{
  \item GPD (PickandsEstimator, medkMAD, medSn, medQn)
  \item GEV (PickandsEstimator)
  \item Pareto (Cramér-von-Mises-Minimum-Distance-Estimator)
  \item Weibull (the quantile based estimator of Boudt/Caliskan/Croux)
  \item Gamma (Cramér-von-Mises-Minimum-Distance-Estimator)
}

For all these families, of course, MLEs and Minimum-Distance-Estimators
are also available through package "distrMod". %\cr
}
\section{Diagnostics}{

We bridge to the diagnostics provided by package "ismev", i.e. our
return objects can be plugged into the diagnostics of this package.\cr

We have the usual diagnostic plots from package "RobAStBase", i.e.
\itemize{
  \item Outylingness plots  \code{outlyingPlotIC}
  \item IC plots \code{plot}
  \item Information plots via \code{infoPlot}
  \item IC comparison plots via \code{comparePlot}
  \item Cniperpoint plots (from package "ROptEst") via \code{CniperPointPlot}
}
but also (adopted from package "distrMod")
\itemize{
  \item qqplots (with confidence bands) via \code{qqplot}
  \item returnlevel plots via \code{returnlevelplot}
}%\cr
}

\section{Starting Point}{
 As a starting point you may look at the included script
  \file{"RobFitsAtRealData.R"} in the scripts folder of the package,
  accessible by
    \code{file.path(system.file(package="RobExtremes"),
             "scripts/RobFitsAtRealData.R")}.
%\cr
}

\details{
\tabular{ll}{
Package: \tab RobExtremes \cr
Version: \tab 1.1.0 \cr
Date: \tab 2018-08-03 \cr
Title: \tab Optimally Robust Estimation for Extreme Value Distributions\cr
Description: \tab Optimally robust estimation for extreme value distributions
using S4 classes and methods \cr
\tab (based on packages distr, distrEx, distrMod, RobAStBase, and ROptEst). \cr
Depends:\tab R (>= 2.14.0), methods, distrMod(>= 2.5.2), ROptEst(>= 1.0), robustbase(>= 0.8-0),
        evd \cr
Suggests: \tab RUnit (>= 0.4.26), ismev (>= 1.39) \cr
Imports: \tab actuar, RobAStRDA, distr, distrEx, RandVar, RobAStBase, startupmsg \cr
Authors:
\tab Bernhard Spangl [contributed smoothed grid values of the Lagrange multipliers]\cr
\tab Sascha Desmettre [contributed smoothed grid values of the Lagrange multipliers]\cr
\tab Eugen Massini [contributed an interactive smoothing routine for smoothing the\cr
\tab       Lagrange multipliers and smoothed grid values of the Lagrange multipliers] \cr
\tab Daria Pupashenko [contributed MDE-estimation for GEV distribution in
the framework of\cr
\tab       her PhD thesis 2011--14]\cr
\tab Gerald Kroisandt [contributed testing routines]\cr
\tab Nataliya Horbenko ["aut","cph"] \cr
\tab Matthias Kohl ["aut", "cph"]\cr
\tab Peter Ruckdeschel ["cre", "aut", "cph"],\cr
Contact: \tab peter.ruckdeschel@uni-oldenburg.de\cr
ByteCompile: \tab yes \cr
LazyLoad: \tab yes \cr
License: \tab LGPL-3 \cr
URL: \tab http://robast.r-forge.r-project.org/\cr
Encoding: \tab latin1 \cr
VCS/SVNRevision: \tab 1091 \cr
}
}

\section{Classes}{
\preformatted{

[*]: there is a generating function with the same name in RobExtremes
[**]:  generating function from distrMod, but with (speeded-up)
       opt.rob-estimators in RobExtremes

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
|>|>|>|>|>"GammaFamily"                 [**]
|>|>|>|>|>"GParetoFamily"               [*]
|>|>|>|>|>"GEVFamily"                   [*]
|>|>|>|>|>"WeibullFamily"               [**]
|>|>|>|>"L2LocationScaleUnion"  /VIRTUAL/     (from distrMod)
|>|>|>|>|>"L2LocationFamily"                  (from distrMod)
|>|>|>|>|>|>"GumbelLocationFamily"      [*]
|>|>|>|>"L2LocScaleShapeUnion"  /VIRTUAL/     (from distrMod)
|>|>|>|>|>"GEVFamilyMuUnknown"          [*]
}
}

\section{Functions}{

\preformatted{
LDEstimator     Estimators for scale-shape models based on
                location and dispersion
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
var                     Generic functions for the computation of functionals
IQR                     Generic functions for the computation of functionals
median                  Generic functions for the computation of functionals
skewness                Generic functions for the computation of functionals
kurtosis                Generic functions for the computation of functionals
Sn                      Generic function for the computation of (conditional)
                        expectations
Qn                      Generic functions for the computation of functionals

}
}
\section{Constants}{
\preformatted{

EULERMASCHERONICONSTANT
APERYCONSTANT

}}

\section{Acknowledgement}{
This package is joint work by Peter Ruckdeschel, Matthias Kohl, and
Nataliya Horbenko (whose PhD thesis went into this package to a large extent),
with contributions by Dasha Pupashenko, Misha Pupashenko, Gerald Kroisandt,
Eugen Massini, Sascha Desmettre, and Bernhard Spangl, in the framework of
project "Robust Risk Estimation" (2011-2016) funded by Volkswagen foundation
(and gratefully ackknowledged). Thanks also goes to the maintainers of CRAN,
in particully to Uwe Ligges who greatly helped us with finding an appropriate
way to store the database of interpolating functions which allow the speed up
-- this is now package RobAStRDA on CRAN. %\cr
}

\author{
Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de},\cr
Matthias Kohl \email{Matthias.Kohl@stamats.de}, and \cr
Nataliya Horbenko \email{nhorbenko@gmail.com},\cr

\emph{Maintainer:}  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
}
\references{
M. Kohl (2005): \emph{Numerical Contributions to the Asymptotic
Theory of Robustness.} PhD Thesis. Bayreuth. Available as
\url{http://r-kurs.de/RRlong.pdf}

P. Ruckdeschel, M. Kohl, T. Stabla, F. Camphausen (2006):
S4 Classes for Distributions, \emph{R News}, \emph{6}(2), 2-6. 
\url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-2.pdf}

M. Kohl, P. Ruckdeschel, H. Rieder (2010):
Infinitesimally Robust Estimation in General Smoothly Parametrized Models.
\emph{Stat. Methods Appl.}, \bold{19}, 333--354.\cr

Ruckdeschel, P. and Horbenko, N. (2011): Optimally-Robust Estimators in Generalized
Pareto Models. \emph{Statistics}. \bold{47}(4), 762--791.\cr

Ruckdeschel, P. and Horbenko, N. (2012): Yet another breakdown point notion:
EFSBP --illustrated at scale-shape models. \emph{Metrika}, \bold{75}(8),
1025--1047.

%a more detailed manual for \pkg{distr}, \pkg{distrSim}, \pkg{distrTEst}, and \pkg{RobExtremes} may be downloaded from
%\url{http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf}\cr

A vignette for packages \pkg{distr}, \pkg{distrSim}, \pkg{distrTEst},
and \pkg{RobExtremes} is included into the mere documentation package \pkg{distrDoc}
and may be called by \code{require("distrDoc");vignette("distr")}.

A homepage to this package is available under \url{http://robast.r-forge.r-project.org/}.

}

\section{Start-up-Banner}{
You may suppress the start-up banner/message completely by setting
\code{options("StartupBanner"="off")} somewhere before loading this package by
\code{library} or \code{require} in your R-code / R-session.
%
If option \code{"StartupBanner"} is not defined (default) or setting
\code{options("StartupBanner"=NULL)} or
\code{options("StartupBanner"="complete")} the complete start-up banner is
displayed.
%
For any other value of option \code{"StartupBanner"} (i.e., not in
\code{c(NULL,"off","complete")}) only the version information is displayed.
%
The same can be achieved by wrapping the \code{library} or \code{require}  call
into either \code{suppressStartupMessages()} or
\code{onlytypeStartupMessages(.,atypes="version")}.
%
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
