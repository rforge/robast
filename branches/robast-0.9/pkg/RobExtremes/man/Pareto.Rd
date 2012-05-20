\name{Pareto}
\alias{Pareto}

\title{Generating function for Pareto-class}
\description{
  Generates an object of class \code{"Pareto"}.
}
\usage{Pareto(shape = 1, Min = 1)}
\arguments{
  \item{shape}{ positive real number: shape parameter of 
    the Pareto distribution. }
  \item{Min}{ positive real number: Min parameter
    of the Pareto distribution }
}
%\details{}
\value{Object of class \code{"Pareto"}}
%\references{}
\author{Nataliya Horbenko \email{Nataliya.Horbenko@itwm.fraunhofer.de}}
\note{The class \code{"Pareto"} is based on the code provided 
  by the package \pkg{actuar} by  Vincent Goulet and Mathieu Pigeon.}
\seealso{\code{\link{Pareto-class}}, \code{\link[actuar:SingleParameterPareto]{dpareto1}}}
\examples{
(P1 <- Pareto(shape = 1, Min = 1))
plot(P1)

E(Pareto()) 
E(P1, function(x){x^2})

## The function is currently defined as
function(shape = 1, Min = 1) 
               new("Pareto", shape = shape, Min = Min)
}

\concept{Pareto}
\keyword{distribution}
\concept{absolutely continuous distribution}
\concept{Pareto distribution}
\concept{generating function}
