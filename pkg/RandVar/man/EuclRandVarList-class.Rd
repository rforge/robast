\name{EuclRandVarList-class}
\docType{class}
\alias{EuclRandVarList-class}
\alias{numberOfMaps}
\alias{numberOfMaps,EuclRandVarList-method}
\alias{dimension,EuclRandVarList-method}
\alias{evalRandVar,EuclRandVarList,numeric,missing-method}
\alias{evalRandVar,EuclRandVarList,matrix,missing-method}
\alias{evalRandVar,EuclRandVarList,numeric,Distribution-method}
\alias{evalRandVar,EuclRandVarList,matrix,Distribution-method}
\alias{imageDistr,EuclRandVarList,Distribution-method}
\alias{t,EuclRandVarList-method}
\alias{show,EuclRandVarList-method}
\alias{Arith,numeric,EuclRandVarList-method}
\alias{Arith,EuclRandVarList,numeric-method}
\alias{Arith,EuclRandVarList,EuclRandVarList-method}
\alias{Math,EuclRandVarList-method}
\alias{\%m\%}
\alias{\%m\%,EuclRandVarList,EuclRandVarList-method}
\alias{\%*\%,matrix,EuclRandVarList-method}
\alias{\%*\%,EuclRandVarList,matrix-method}
\alias{E,UnivariateDistribution,EuclRandVarList,missing-method}
\alias{E,AbscontDistribution,EuclRandVarList,missing-method}
\alias{E,DiscreteDistribution,EuclRandVarList,missing-method}
\alias{E,MultivariateDistribution,EuclRandVarList,missing-method}
\alias{E,DiscreteMVDistribution,EuclRandVarList,missing-method}
\alias{E,UnivariateCondDistribution,EuclRandVarList,numeric-method}
\alias{E,AbscontCondDistribution,EuclRandVarList,numeric-method}
\alias{E,DiscreteCondDistribution,EuclRandVarList,numeric-method}

\title{List of Euclidean random variables}
\description{Create a list of Euclidean random variables}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("EuclRandVarList", ...)}.
  More frequently they are created via the generating function 
  \code{EuclRandVarList}. 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}. A list of Euclidean
      random variables. }
  }
}
\section{Extends}{
Class \code{"list"}, from data part.\cr
Class \code{"vector"}, by class \code{"list"}.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "EuclRandVariable", to = "EuclRandVarList")}: 
      create a \code{"EuclRandVarList"} object from a Euclidean random variable. }

    \item{coerce}{\code{signature(from = "EuclRandMatrix", to = "EuclRandVarList")}: 
      create a  \code{"EuclRandVarList"} object from a Euclidean random matrix. }

    \item{numberOfMaps}{\code{signature(object = "EuclRandVarList")}:
      number of functions contained in the slots \code{Map} of
      the members of \code{object}. }

    \item{dimension}{\code{signature(object = "EuclRandVarList")}: 
      dimension of the Euclidean random variable. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandVarList", x = "numeric")}:
      evaluate the elements of \code{RandVar} at \code{x}. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandVarList", x = "matrix")}:
      evaluate the elements of \code{RandVar} at rows of \code{x}. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandVarList", x = "numeric", distr = "Distribution")}:
      evaluate the elements of \code{RandVar} at \code{x} assuming 
      a probability space with distribution \code{distr}. In case \code{x}
      does not lie in the support of \code{distr} \code{NA} is returned. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandVarList", x = "matrix", distr = "Distribution")}:
      evaluate the elements of \code{RandVar} at rows of \code{x}
      assuming a probability space with distribution \code{distr}. For those
      rows of \code{x} which do not lie in the support of \code{distr} 
      \code{NA} is returned. }

    \item{imageDistr}{\code{signature(RandVar = "EuclRandVarList", distr = "Distribution")}: 
      image distribution of \code{distr} under \code{RandVar}. Returns
      an object of class \code{"DistrList"}. }

    \item{show}{\code{signature(object = "EuclRandVarList")} }

    \item{t}{\code{signature(x = "EuclRandVarList")}: 
      returns an object of class \code{"EuclRandVarList"} where the
      rhe results of the functions in the slots \code{Map} of the members of 
      \code{x} are transposed. }

    \item{\%m\%}{\code{signature(x = "EuclRandVarList", y = "EuclRandVarList")}: 
      matrix multiplication for objects of class \code{"EuclRandVarList"}. 
     Generates an object of class  \code{"EuclRandVarList"}.}

    \item{\%*\%}{\code{signature(x = "matrix", y = "EuclRandVarList")}:  
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{\%*\%}{\code{signature(x = "EuclRandVarList", y = "matrix")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{Arith}{\code{signature(e1 = "numeric", e2 = "EuclRandVarList")}:
      Given a numeric vector \code{e1}, a list of Euclidean random variables \code{e2} 
      and an arithmetic operator \code{op}, the list of Euclidean random variables 
      \code{e1 op e2} is returned. }

    \item{Arith}{\code{signature(e1 = "EuclRandVarList", e2 = "numeric")}:
      Given a numeric vector \code{e2}, a list of Euclidean random variables \code{e1} 
      and an arithmetic operator \code{op}, the list of Euclidean random variables 
      \code{e1 op e2} is returned. }

    \item{Arith}{\code{signature(e1 = "EuclRandVarList", e2 = "EuclRandVarList")}:
      Given two lists of Euclidean random variables \code{e1}, \code{e2} and an 
      arithmetic operator \code{op}, the list of Euclidean random variables 
      \code{e1 op e2} is returned. }

    \item{Math}{\code{signature(x = "EuclRandVarList")}:
      Given a \code{"Math"} group generic \code{fct}, the list of Euclidean random 
      variables \code{fct(x)} is returned. }

    \item{E}{\code{signature(object = "UnivariateDistribution", fun = "EuclRandVarList", cond = "missing")}: 
      expectation of \code{fun} under univariate distributions. }

    \item{E}{\code{signature(object = "AbscontDistribution", fun = "EuclRandVarList", cond = "missing")}:
      expectation of \code{fun} under absolutely continuous univariate distributions. }

    \item{E}{\code{signature(object = "DiscreteDistribution", fun = "EuclRandVarList", cond = "missing")}:
      expectation of \code{fun} under discrete univariate distributions. }

    \item{E}{\code{signature(object = "MultivariateDistribution", fun = "EuclRandVarList", cond = "missing")}: 
      expectation of \code{fun} under multivariate distributions. }

    \item{E}{\code{signature(object = "DiscreteMVDistribution", fun = "EuclRandVarList", cond = "missing")}: 
      expectation of \code{fun} under discrete multivariate distributions. }

    \item{E}{\code{signature(object = "UnivariateCondDistribution", fun = "EuclRandVarList", cond = "numeric")}: 
      expectation of \code{fun} under conditional univariate distributions. }

    \item{E}{\code{signature(object = "AbscontCondDistribution", fun = "EuclRandVarList", cond = "numeric")}: 
      expectation of \code{fun} under absolutely continuous conditional univariate distributions. }

    \item{E}{\code{signature(object = "DiscreteCondDistribution", fun = "EuclRandVarList", cond = "numeric")}: 
      expectation of \code{fun} under discrete conditional univariate distributions. }
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{EuclRandMatrix}}, \code{\link{RandVariable-class}}, 
    \code{\link{EuclRandVariable-class}}, \code{\link{EuclRandMatrix-class}},
    \code{\link[distr]{Distribution-class}}, \code{\link[base:groupGeneric]{Arith}}, 
    \code{\link[base:groupGeneric]{Math}}, \code{\link[distrEx]{E}}}
\examples{
L1 <- list(function(x){x}, function(x){x^2}, function(x){x^3}, function(x){x^4}, 
           function(x){x^5}, function(x){x^6})
L2 <- list(function(x){exp(x)}, function(x){abs(x)}, 
           function(x){sin(x)}, function(x){floor(x)})

R1 <- new("EuclRandVariable", Map = L2, Domain = Reals(), Range = Reals())
R2 <- EuclRandMatrix(Map = L1, ncol = 2, Domain = Reals(), dimension = 1)
R3 <- EuclRandMatrix(Map = L2, ncol = 2, Domain = Reals(), dimension = 1)

(RL1 <- new("EuclRandVarList", list(R1, R2, R3)))
dimension(RL1)
as(R1, "EuclRandVarList")
as(R2, "EuclRandVarList")

Map(exp(RL1)[[1]]) # "Math" group

## "Arith" group
Map((1 + RL1)[[1]])
Map((RL1 * 2)[[2]])
Map((RL1 / RL1)[[3]])
}
\concept{random variable}
\concept{random matrix}
\keyword{classes}
