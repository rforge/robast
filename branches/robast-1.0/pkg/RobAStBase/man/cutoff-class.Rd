\name{cutoff-class}
\docType{class}
\alias{cutoff-class}
\alias{cutoff.quantile<-,cutoff-method}
\alias{cutoff.quantile<-}
\alias{cutoff.quantile,cutoff-method}
\alias{cutoff.quantile}
\alias{name,cutoff-method}
\alias{fct,cutoff-method}

\title{Cutoff class for distance-distance plots}
\description{Class of methods to determine cutoff point for distance-distance plots;
  used to derive other cutoff methods later by method dispatch.
}
\section{Objects from the Class}{
  Objects could in principle be created by calls of the form \code{new("cutoff", ...)}.
  More frequently they are created via the generating function 
  \code{\link{cutoff}}, respectively via the helper functions 
  \code{cutoff.sememp} and \code{cutoff.chisq}. 
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{ object of class \code{"character"}; defaults to \code{"empirical"}
                         in prototype; }

    \item{\code{fct}:}{ an object of of class \code{"function"}; 
                        for this class layer, this function 
                         must only have one argument \code{data} (which may
    but need not be used to determine the cutoff point empirically);
    in derived classes this restriction could be dropped, if corresponding
    special methods for \code{ddPlot} are derived. Defaults to 
    \code{function(data) quantile(data)}.}

    \item{\code{cutoff.quantile}:}{ Object of class \code{"numeric"}:
      a probability (in [0,1]) to determine the respective quantile 
      (empirical or theoretical) to plot the cutoff line; defaults
      to \code{0.95} in prototype;}

  }
}
\section{Methods}{
  \describe{

    \item{cutoff.quantile}{\code{signature(object = "cutoff")}: 
      accessor function for slot \code{cutoff.quantile}. }

    \item{cutoff.quantile<-}{\code{signature(object = "cutoff")}: 
      replacement function for slot \code{cutoff.quantile}. }

    \item{fct}{\code{signature(object = "cutoff")}: 
      accessor function for slot \code{fct}. }

    \item{name}{\code{signature(object = "cutoff")}: 
      accessor function for slot \code{name}. }


   }
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{ddPlot}}, \code{\link{outlyingPlotIC}} \code{\link{cutoff}}}
\examples{
cutoff()
}
\keyword{classes}
