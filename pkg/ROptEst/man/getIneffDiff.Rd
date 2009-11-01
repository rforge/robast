\name{getIneffDiff}
\alias{getIneffDiff}
\alias{getIneffDiff-methods}
\alias{getIneffDiff,numeric,L2ParamFamily,UncondNeighborhood,asMSE-method}

\title{Generic Function for the Computation of Inefficiency Differences}
\description{
  Generic function for the computation of inefficiency differencies.
  This function is rarely called directly. It is used to compute
  the radius minimax IC and the least favorable radius.
}
\usage{
getIneffDiff(radius, L2Fam, neighbor, risk, ...)

\S4method{getIneffDiff}{numeric,L2ParamFamily,UncondNeighborhood,asMSE}(
          radius, L2Fam, neighbor, risk, loRad, upRad, loRisk, upRisk, 
          z.start = NULL, A.start = NULL, upper.b = NULL, lower.b = NULL,
          OptOrIter = "iterate", MaxIter, eps, warn, loNorm = NULL, upNorm = NULL,
          verbose = NULL, ...)
}
\arguments{
  \item{radius}{ neighborhood radius. }
  \item{L2Fam}{ L2-differentiable family of probability measures. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{loRad}{ the lower end point of the interval to be searched. }
  \item{upRad}{ the upper end point of the interval to be searched. }
  \item{loRisk}{ the risk at the lower end point of the interval. }
  \item{upRisk}{ the risk at the upper end point of the interval. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{upper.b}{ upper bound for the optimal clipping bound. }
  \item{lower.b}{ lower bound for the optimal clipping bound. }
  \item{OptOrIter}{character; which method to be used for determining Lagrange
  multipliers \code{A} and \code{a}: if (partially) matched to \code{"optimize"},
  \code{getLagrangeMultByOptim} is used; otherwise: by default, or if matched to
  \code{"iterate"} or to \code{"doubleiterate"},
  \code{getLagrangeMultByIter} is used. More specifically,
  when using \code{getLagrangeMultByIter}, and if argument \code{risk} is of
  class \code{"asGRisk"}, by default and if matched to \code{"iterate"}
  we use only one (inner) iteration, if matched to \code{"doubleiterate"}
  we use up to \code{Maxiter} (inner) iterations.}
  \item{MaxIter}{ the maximum number of iterations }
  \item{eps}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }  
  \item{loNorm}{object of class \code{"NormType"}; used in selfstandardization
                to evaluate the bias of the current IC in the norm of the lower 
                bound}  
  \item{upNorm}{object of class \code{"NormType"}; used in selfstandardization
                to evaluate the bias of the current IC in the norm of the upper 
                bound}  
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed }
  \item{\dots}{further arguments to be passed on to \code{getInfRobIC}}
}
%\details{}
\value{The inefficieny difference between the left and
  the right margin of a given radius interval is computed.
}
\section{Methods}{
\describe{
  \item{radius = "numeric", L2Fam = "L2ParamFamily", 
        neighbor = "UncondNeighborhood", risk = "asMSE":}{ 
     computes difference of asymptotic MSE--inefficiency for
     the boundaries of a given radius interval.}
}}
\references{ 
  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications, \emph{17}(1) 13-40.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2001) The Costs of not Knowing
  the Radius. Submitted. Appeared as discussion paper Nr. 81. 
  SFB 373 (Quantification and Simulation of Economic Processes),
  Humboldt University, Berlin; also available under
  \url{www.uni-bayreuth.de/departments/math/org/mathe7/RIEDER/pubs/RR.pdf}

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{radiusMinimaxIC}}, \code{\link{leastFavorableRadius}}}
%\examples{}
\keyword{robust}
