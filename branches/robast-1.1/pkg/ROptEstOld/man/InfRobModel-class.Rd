\name{InfRobModel-class}
\docType{class}
\alias{InfRobModel-class}
\alias{neighbor<-,InfRobModel-method}
\alias{show,InfRobModel-method}

\title{Robust model with infinitesimal (unconditional) neighborhood}
\description{Class of robust models with infinitesimal (unconditional) neighborhoods;
  i.e., the neighborhood is shrinking at a rate of \eqn{\sqrt{n}}{sqrt(n)}.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("InfRobModel", ...)}.
  More frequently they are created via the generating function 
  \code{InfRobModel}. 
}
\section{Slots}{
  \describe{
    \item{\code{center}:}{Object of class \code{"ProbFamily"}. }
    \item{\code{neighbor}:}{Object of class \code{"UncondNeighborhood"}. }
  }
}
\section{Extends}{
Class \code{"RobModel"}, directly.
}
\section{Methods}{
  \describe{
    \item{neighbor<-}{\code{signature(object = "InfRobModel")}: 
      replacement function for slot \code{neighbor<-} }

    \item{show}{\code{signature(object = "InfRobModel")}}
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ProbFamily-class}}, \code{\link{UncondNeighborhood-class}},
  \code{\link{InfRobModel}}}
\examples{
new("InfRobModel")
}
\concept{robust model}
\concept{infinitesimal neighborhood}
\concept{shrinking neighborhood}
\keyword{classes}
\keyword{models}
