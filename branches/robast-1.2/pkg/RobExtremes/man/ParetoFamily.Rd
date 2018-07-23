\name{ParetoFamily}
\alias{ParetoFamily}

\title{Generating function for Generalized Pareto families}
\description{
  Generates an object of class \code{"ParetoFamily"} which
  represents a Pareto family.
}
\usage{
ParetoFamily(Min = 1, shape = 0.5, trafo = NULL, start0Est = NULL,
                    withCentL2 = FALSE)
}
\arguments{
  \item{Min}{ real: known/fixed threshold/location parameter }
  \item{shape}{ positive real: shape parameter }
  \item{trafo}{ matrix or NULL: transformation of the parameter }
  \item{start0Est}{ startEstimator --- if \code{NULL}
          \code{log(2)/log(median/Min)} is used }
  \item{withCentL2}{logical: shall L2 derivative be centered by substracting
       the E()? Defaults to \code{FALSE}, but higher accuracy can be achieved
       when set to \code{TRUE}.}
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
}
\value{Object of class \code{"ParetoFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.\cr

M.~Kohl, P. Ruckdeschel, H.~Rieder (2010):
Infinitesimally Robust Estimation in General Smoothly Parametrized Models.
\emph{Stat. Methods Appl.}, \bold{19}, 333--354.\cr

Ruckdeschel, P. and Horbenko, N. (2011): Optimally-Robust Estimators in Generalized
Pareto Models. ArXiv 1005.1476. To appear at \emph{Statistics}.
DOI: 10.1080/02331888.2011.628022. \cr

Ruckdeschel, P. and Horbenko, N. (2012): Yet another breakdown point notion:
EFSBP --illustrated at scale-shape models. \emph{Metrika}, \bold{75}(8),
1025--1047.
  }

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}\cr
        Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}\cr
        Nataliya Horbenko \email{nhorbenko@gmail.com}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\linkS4class{Pareto}}}
\examples{
(P1 <- ParetoFamily())
FisherInfo(P1)
checkL2deriv(P1)
}
\concept{Pareto model}
\keyword{models}
