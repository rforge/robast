\name{GumbelLocationFamily}
\alias{GumbelLocationFamily}

\title{Generating function for Gumbel location families}
\description{
  Generates an object of class \code{"L2LocationFamily"} which
  represents a Gumbel location family.
}
\usage{
GumbelLocationFamily(loc = 0, scale = 1, trafo)
}
\arguments{
  \item{loc}{ location parameter }
  \item{scale}{ scale parameter }
  \item{trafo}{ function in \code{param} or matrix: transformation of the parameter }
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
}
\value{Object of class \code{"L2LocationFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link{Gumbel-class}}}
\examples{
##current implementation is:
theta <- 0
names(theta) <- "loc"
GL <- ParamFamily(name = "Gumbel location family",
          param = ParamFamParameter(name = "location parameter", main = theta),
          startPar = function(x,...) c(min(x),max(x)),
          distribution = Gumbel(loc = 0, scale = 1), ## scale known!
          modifyParam = function(theta){ Gumbel(loc = theta, scale = 1) },
          props = paste(c("The Gumbel location family is invariant under",
                    "the group of transformations 'g(x) = x + loc'",
                    "with location parameter 'loc'"), collapse = " "))
GL

(G1 <- GumbelLocationFamily())
plot(G1)
Map(L2deriv(G1)[[1]])
checkL2deriv(G1)
}
\concept{Gumbel location model}
\concept{location model}
\keyword{models}
