\name{finiteSampleCorrection}
\Rdversion{1.1}
\alias{finiteSampleCorrection}
\title{Function to compute finite-sample corrected radii}
\description{
Given some radius and some sample size the function computes
the corresponding finite-sample corrected radius.
}
\usage{
finiteSampleCorrection(r, n, model = "locsc")
}
\arguments{
  \item{r}{ asymptotic radius (non-negative numeric) }
  \item{n}{ sample size }
  \item{model}{ has to be \code{"locsc"} (for location and scale),
  \code{"loc"} (for location) or \code{"sc"} (for scale), respectively. }
}
\details{
The finite-sample correction is based on empirical results obtained via 
simulation studies. 

Given some radius of a shrinking contamination neighborhood which leads 
to an asymptotically optimal robust estimator, the finite-sample empirical 
MSE based on contaminated samples was minimized for this class of 
asymptotically optimal estimators and the corresponding finite-sample 
radius determined and saved.

The computation is based on the saved results of these Monte-Carlo simulations.
}
\value{Finite-sample corrected radius.}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.
  Extended version: \url{http://r-kurs.de/RRlong.pdf}
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{roblox}}, \code{\link{rowRoblox}}, 
         \code{\link{colRoblox}} }
\examples{
finiteSampleCorrection(n = 3, r = 0.001, model = "locsc")
finiteSampleCorrection(n = 10, r = 0.02, model = "loc")
finiteSampleCorrection(n = 250, r = 0.15, model = "sc")
}
\concept{normal location}
\concept{normal scale}
\concept{normal location and scale}
\concept{finite-sample correction}
\keyword{robust}
