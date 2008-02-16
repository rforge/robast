\name{lowerCaseRadius}
\alias{lowerCaseRadius}
\alias{lowerCaseRadius-methods}
\alias{lowerCaseRadius,L2ParamFamily,ContNeighborhood,asMSE,BiasType-method}
\alias{lowerCaseRadius,L2ParamFamily,TotalVarNeighborhood,asMSE,BiasType-method}

\title{Computation of the lower case radius}
\description{
  The lower case radius is computed; confer Subsection 2.1.2 
  in Kohl (2005) and formula (4.5) in Ruckdeschel (2005).
}
\usage{
lowerCaseRadius(L2Fam, neighbor, risk, biastype, ...)
}
\arguments{
  \item{L2Fam}{ L2 differentiable parametric family }
  \item{neighbor}{ object of class \code{"Neighborhood"} }
  \item{risk}{ object of class \code{"RiskType"} }
  \item{biastype}{ object of class \code{"BiasType"} }
  \item{\dots}{ additional parameters }
}
%\details{}
\value{lower case radius}
\section{Methods}{
\describe{
  \item{L2Fam = "L2ParamFamily", neighbor = "ContNeighborhood", risk = "asMSE",
        biastype = "BiasType"}{
    lower case radius for risk \code{"asMSE"} in case of \code{"ContNeighborhood"}
    for symmetric bias.}

  \item{L2Fam = "L2ParamFamily", neighbor = "TotalVarNeighborhood", risk = "asMSE",
        biastype = "BiasType"}{
    lower case radius for risk \code{"asMSE"} in case of \code{"TotalVarNeighborhood"};
    (argument biastype is just for signature reasons).}
}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},
        Peter Ruckdeschel \email{Peter.Ruckdeschel@uni-bayreuth.de}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\link[RobAStBase]{Neighborhood-class}}}
\examples{
lowerCaseRadius(BinomFamily(size = 10), ContNeighborhood(), asMSE())
lowerCaseRadius(BinomFamily(size = 10), TotalVarNeighborhood(), asMSE())
}
\concept{lower case radius}
\keyword{}
