\name{.getImageDistr}
\alias{.getImageDistr}

\title{Help function for generation of image distributions}
\description{
  Generates an imageDistribution f(distr).
}
\usage{
.getImageDistr(f, distr)
}
\arguments{
  \item{f}{a function with values in R}
  \item{distr}{an object of class \code{"Abscontdistribution"} }
}
\value{image distribution of \code{distr} under \code{f}}
\details{works even if \code{f(distr)} fails; 
  if anything else fails does simulations of \code{f(x)}, 
  \code{x} according to \code{distr}; uses \code{RtoDPQ} then;
  does not check whether \code{f(distr)} has point masses}
%\references{}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
%\examples{}
\concept{image distribution}
\concept{random variable}
\keyword{internal}
