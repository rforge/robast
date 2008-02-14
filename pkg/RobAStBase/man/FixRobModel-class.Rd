\name{FixRobModel-class}
\docType{class}
\alias{FixRobModel-class}
\alias{neighbor<-,FixRobModel-method}
\alias{show,FixRobModel-method}

\title{Robust model with fixed (unconditional) neighborhood}
\description{Class of robust models with fixed (unconditional) neighborhoods.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("FixRobModel", ...)}.
  More frequently they are created via the generating function 
  \code{FixRobModel}. 
}
\section{Slots}{
  \describe{
    \item{\code{center}:}{Object of class \code{"ProbFamily"}. }
    \item{\code{neighbor}:}{Object of class \code{"UncondNeighborhood"}.}
  }
}
\section{Extends}{
Class \code{"RobModel"}, directly.
}
\section{Methods}{
  \describe{
    \item{neighbor<-}{\code{signature(object = "FixRobModel")}: 
      replacement function for slot \code{neighbor<-} }

    \item{show}{\code{signature(object = "FixRobModel")}}
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
  \code{\link{FixRobModel}}}
\examples{
new("FixRobModel")
}
\concept{robust model}
\concept{fixed neighborhood}
\keyword{classes}
\keyword{models}
