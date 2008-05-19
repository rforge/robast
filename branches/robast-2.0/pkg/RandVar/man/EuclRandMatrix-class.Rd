\name{EuclRandMatrix-class}
\docType{class}
\alias{EuclRandMatrix-class}
\alias{coerce,EuclRandMatrix,EuclRandVarList-method}
\alias{Dim}
\alias{Dim,EuclRandMatrix-method}
\alias{Dim<-}
\alias{Dim<-,EuclRandMatrix-method}
\alias{[,EuclRandMatrix-method}
\alias{ncol,EuclRandMatrix-method}
\alias{nrow,EuclRandMatrix-method}
\alias{dimension,EuclRandMatrix-method}
\alias{evalRandVar,EuclRandMatrix,numeric,missing-method}
\alias{evalRandVar,EuclRandMatrix,matrix,missing-method}
\alias{evalRandVar,EuclRandMatrix,numeric,Distribution-method}
\alias{evalRandVar,EuclRandMatrix,matrix,Distribution-method}
\alias{t,EuclRandMatrix-method}
\alias{show,EuclRandMatrix-method}
\alias{\%*\%,matrix,EuclRandMatrix-method}
\alias{\%*\%,numeric,EuclRandMatrix-method}
\alias{\%*\%,EuclRandMatrix,matrix-method}
\alias{\%*\%,EuclRandMatrix,numeric-method}
\alias{\%*\%,EuclRandMatrix,EuclRandMatrix-method}
\alias{Arith,numeric,EuclRandMatrix-method}
\alias{Arith,EuclRandMatrix,numeric-method}
\alias{Arith,EuclRandMatrix,EuclRandMatrix-method}
\alias{Math,EuclRandMatrix-method}
\alias{E,UnivariateDistribution,EuclRandMatrix,missing-method}
\alias{E,AbscontDistribution,EuclRandMatrix,missing-method}
\alias{E,DiscreteDistribution,EuclRandMatrix,missing-method}
\alias{E,MultivariateDistribution,EuclRandMatrix,missing-method}
\alias{E,DiscreteMVDistribution,EuclRandMatrix,missing-method}
\alias{E,UnivariateCondDistribution,EuclRandMatrix,numeric-method}
\alias{E,AbscontCondDistribution,EuclRandMatrix,numeric-method}
\alias{E,DiscreteCondDistribution,EuclRandMatrix,numeric-method}

\title{Euclidean random matrix}
\description{Class of Euclidean random matrices.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("EuclRandMatrix", ...)}.
  More frequently they are created via the generating function 
  \code{EuclRandMatrix}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}:}{ vector of positive integers: Dimensions of the random matrix. }
    \item{\code{Map}:}{ Object of class \code{"list"}: list of functions. }
    \item{\code{Domain}:}{ Object of class \code{"OptionalrSpace"} 
      domain of the random matrix. }
    \item{\code{Range}:}{ Object of class \code{"OptionalrSpace"} 
      range of the random matrix. }
  }
}
\section{Extends}{
Class \code{"EuclRandVariable"}, directly.\cr
Class \code{"RandVariable"}, by class \code{"EuclRandVariable"}.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "EuclRandMatrix", to = "EuclRandVarList")}: 
      create a \code{"EuclRandVarList"} object from a Euclidean random matrix. }

    \item{[}{\code{signature(x = "EuclRandMatrix")}: generates
      a new Euclidean random variable/matrix by extracting elements of
      the slot \code{Map} of \code{x}. }

    \item{Dim}{\code{signature(object = "EuclRandMatrix")}: accessor function
      for slot \code{Dim}. }

    \item{Dim<-}{\code{signature(object = "EuclRandMatrix", )}: replacement
      function for slot \code{Dim}. }

    \item{ncol}{\code{signature(x = "EuclRandMatrix")}: number of columns of \code{x}. }

    \item{nrow}{\code{signature(x = "EuclRandMatrix")}: number of rows of \code{x}. }

    \item{dimension}{\code{signature(object = "EuclRandMatrix")}: dimension 
      of the Euclidean random variable. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandMatrix", x = "numeric")}: 
      evaluate the slot \code{Map} of \code{RandVar} at \code{x}. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandMatrix", x = "matrix")}:
      evaluate the slot \code{Map} of \code{RandVar} at \code{x}. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandMatrix", x = "numeric", distr = "Distribution")}:
      evaluate the slot \code{Map} of \code{RandVar} at \code{x} assuming 
      a probability space with distribution \code{distr}. In case \code{x}
      does not lie in the support of \code{distr} \code{NA} is returned. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandMatrix", x = "matrix", distr = "Distribution")}:
      evaluate the slot \code{Map} of \code{RandVar} at rows of \code{x}
      assuming a probability space with distribution \code{distr}. For those
      rows of \code{x} which do not lie in the support of \code{distr} 
      \code{NA} is returned. }

    \item{t}{\code{signature(x = "EuclRandMatrix")}: transposes \code{x}. In
      addition, the results of the functions in the slot \code{Map} of  
      \code{x} are transposed. }

    \item{show}{\code{signature(object = "EuclRandMatrix")} }

    \item{\%*\%}{\code{signature(x = "matrix", y = "EuclRandMatrix")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{\%*\%}{\code{signature(x = "numeric", y = "EuclRandMatrix")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{\%*\%}{\code{signature(x = "EuclRandVariable", y = "EuclRandMatrix")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{\%*\%}{\code{signature(x = "EuclRandMatrix", y = "matrix")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{\%*\%}{\code{signature(x = "EuclRandMatrix", y = "numeric")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{\%*\%}{\code{signature(x = "EuclRandMatrix", y = "EuclRandMatrix")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{\%*\%}{\code{signature(x = "EuclRandMatrix", y = "EuclRandVariable")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{Arith}{\code{signature(e1 = "numeric", e2 = "EuclRandMatrix")}: 
      Given a numeric vector \code{e1}, a Euclidean random matrix \code{e2} 
      and an arithmetic operator \code{op}, the Euclidean random matrix 
      \code{e1 op e2} is returned. }

    \item{Arith}{\code{signature(e1 = "EuclRandMatrix", e2 = "numeric")}: 
      Given a Euclidean random matrix \code{e1}, a numeric vector \code{e2},  
      and an arithmetic operator \code{op}, the Euclidean random matrix 
      \code{e1 op e2} is returned. }

    \item{Arith}{\code{signature(e1 = "EuclRandMatrix", e2 = "EuclRandMatrix")}: 
      Given two Euclidean random matrices \code{e1} and \code{e2},  
      and an arithmetic operator \code{op}, the Euclidean random matrix 
      \code{e1 op e2} is returned. }

    \item{Math}{\code{signature(x = "EuclRandMatrix")}: 
      Given a \code{"Math"} group generic \code{fct}, the Euclidean random 
      matrix \code{fct(x)} is returned. }

    \item{E}{\code{signature(object = "UnivariateDistribution", fun = "EuclRandMatrix", cond = "missing")}: 
      expectation of \code{fun} under univariate distributions. }

    \item{E}{\code{signature(object = "AbscontDistribution", fun = "EuclRandMatrix", cond = "missing")}: 
      expectation of \code{fun} under absolutely continuous univariate distributions. }

    \item{E}{\code{signature(object = "DiscreteDistribution", fun = "EuclRandMatrix", cond = "missing")}: 
      expectation of \code{fun} under discrete univariate distributions. }

    \item{E}{\code{signature(object = "MultivariateDistribution", fun = "EuclRandMatrix", cond = "missing")}: 
      expectation of \code{fun} under multivariate distributions. }

    \item{E}{\code{signature(object = "DiscreteMVDistribution", fun = "EuclRandMatrix", cond = "missing")}: 
      expectation of \code{fun} under discrete multivariate distributions. }

    \item{E}{\code{signature(object = "UnivariateCondDistribution", fun = "EuclRandMatrix", cond = "numeric")}: 
      conditional expectation of \code{fun} under conditional univariate distributions. }

    \item{E}{\code{signature(object = "AbscontCondDistribution", fun = "EuclRandMatrix", cond = "numeric")}: 
      conditional expectation of \code{fun} under absolutely continuous conditional univariate distributions. }

    \item{E}{\code{signature(object = "DiscreteCondDistribution", fun = "EuclRandMatrix", cond = "numeric")}: 
      conditional expectation of \code{fun} under discrete conditional univariate distributions. }
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{EuclRandMatrix}}, \code{\link{RandVariable-class}}, 
    \code{\link{EuclRandVariable-class}}, \code{\link{EuclRandVarList-class}}, 
    \code{\link[distr]{Distribution-class}}, \code{\link[base:groupGeneric]{Arith}}, 
    \code{\link[base:groupGeneric]{Math}}, \code{\link[distrEx]{E}}}
\examples{
L1 <- list(function(x){x}, function(x){x^2}, function(x){x^3}, function(x){x^4}, 
           function(x){x^5}, function(x){x^6})
L2 <- list(function(x){exp(x)}, function(x){abs(x)}, 
           function(x){sin(x)}, function(x){floor(x)})

R1 <- new("EuclRandMatrix", Map = L1, Dim = as.integer(c(3,2)), 
                            Domain = Reals(), Range = Reals())
dimension(R1)
R1[1:2, 2]
R1[1:2, 1:2]
Map(R1[1,2])
Map(t(R1)[2,1])

R2 <- EuclRandMatrix(Map = L2, ncol = 2, Domain = Reals(), dimension = 1)
dimension(R2)
(DL <- imageDistr(R2, Norm()))
plot(DL)

Map(gamma(R2)) # "Math" group

## "Arith" group
Map(2/R1)
Map(R2 * R2)
}
\concept{random variable}
\concept{random matrix}
\keyword{classes}
\keyword{arith}
\keyword{math}
