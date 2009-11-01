\name{lowerCaseRadius}
\alias{lowerCaseRadius}
\alias{lowerCaseRadius-methods}
\alias{lowerCaseRadius,L2ParamFamily,ContNeighborhood,asMSE-method}
\alias{lowerCaseRadius,L2ParamFamily,TotalVarNeighborhood,asMSE-method}

\title{Computation of the lower case radius}
\description{
  The lower case radius is computed; confer Subsection 2.1.2 
  in Kohl (2005).
}
\usage{
lowerCaseRadius(L2Fam, neighbor, risk, ...)
}
\arguments{
  \item{L2Fam}{ L2 differentiable parametric family }
  \item{neighbor}{ object of class \code{"Neighborhood"} }
  \item{risk}{ object of class \code{"RiskType"} }
  \item{\dots}{ additional parameters }
}
%\details{}
\value{lower case radius}
\section{Methods}{
\describe{
  \item{L2Fam = "L2ParamFamily", neighbor = "ContNeighborhood", risk = "asMSE"}{
    lower case radius for risk \code{"asMSE"} in case of \code{"ContNeighborhood"}.}

  \item{L2Fam = "L2ParamFamily", neighbor = "TotalVarNeighborhood", risk = "asMSE"}{
    lower case radius for risk \code{"asMSE"} in case of \code{"TotalVarNeighborhood"}.}
}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link{Neighborhood-class}}}
\examples{
lowerCaseRadius(BinomFamily(size = 10), ContNeighborhood(), asMSE())
lowerCaseRadius(BinomFamily(size = 10), TotalVarNeighborhood(), asMSE())
}
\concept{lower case radius}
\keyword{robust}
