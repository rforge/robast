\name{EuclRandVariable-class}
\docType{class}
\alias{EuclRandVariable-class}
\alias{coerce,EuclRandVariable,EuclRandMatrix-method}
\alias{coerce,EuclRandVariable,EuclRandVarList-method}
\alias{Range<-,EuclRandVariable-method}
\alias{[,EuclRandVariable-method}
\alias{evalRandVar}
\alias{evalRandVar,EuclRandVariable,numeric,missing-method}
\alias{evalRandVar,EuclRandVariable,matrix,missing-method}
\alias{evalRandVar,EuclRandVariable,numeric,Distribution-method}
\alias{evalRandVar,EuclRandVariable,matrix,Distribution-method}
\alias{imageDistr}
\alias{imageDistr,EuclRandVariable,Distribution-method}
\alias{dimension,EuclRandVariable-method}
\alias{t,EuclRandVariable-method}
\alias{\%*\%,matrix,EuclRandVariable-method}
\alias{\%*\%,numeric,EuclRandVariable-method}
\alias{\%*\%,EuclRandVariable,matrix-method}
\alias{\%*\%,EuclRandVariable,numeric-method}
\alias{\%*\%,EuclRandVariable,EuclRandVariable-method}
\alias{\%*\%,EuclRandVariable,EuclRandMatrix-method}
\alias{\%*\%,EuclRandMatrix,EuclRandVariable-method}
\alias{Arith,numeric,EuclRandVariable-method}
\alias{Arith,EuclRandVariable,numeric-method}
\alias{Arith,EuclRandVariable,EuclRandVariable-method}
\alias{Math,EuclRandVariable-method}
\alias{E,UnivariateDistribution,EuclRandVariable,missing-method}
\alias{E,AbscontDistribution,EuclRandVariable,missing-method}
\alias{E,DiscreteDistribution,EuclRandVariable,missing-method}
\alias{E,MultivariateDistribution,EuclRandVariable,missing-method}
\alias{E,DiscreteMVDistribution,EuclRandVariable,missing-method}
\alias{E,UnivariateCondDistribution,EuclRandVariable,numeric-method}
\alias{E,AbscontCondDistribution,EuclRandVariable,numeric-method}
\alias{E,DiscreteCondDistribution,EuclRandVariable,numeric-method}

\title{Euclidean random variable}
\description{Class of Euclidean random variables.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("EuclRandVariable", ...)}.
  More frequently they are created via the generating function 
  \code{EuclRandVariable}. 
}
\section{Slots}{
  \describe{
    \item{\code{Map}}{Object of class \code{"list"}: list of functions. }
    \item{\code{Domain}}{Object of class \code{"OptionalrSpace"}: 
      domain of the random variable. }
    \item{\code{Range}}{Object of class \code{"EuclideanSpace"}:
      range of the random variable. }
  }
}
\section{Extends}{
Class \code{"RandVariable"}, directly.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "EuclRandVariable", to = "EuclRandMatrix")}: 
      create a \code{"EuclRandMatrix"} object from a Euclidean random variable. }

    \item{coerce}{\code{signature(from = "EuclRandVariable", to = "EuclRandVarList")}: 
      create a \code{"EuclRandVarList"} object from a Euclidean random variable. }

    \item{Range<-}{\code{signature(object = "EuclRandVariable")}:
      replacement function for the slot \code{Range}. }

    \item{[}{\code{signature(x = "EuclRandVariable")}: generates
      a new Euclidean random variable by extracting elements of
      the slot \code{Map} of \code{x}.  }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandVariable", x = "numeric", distr = "missing")}:
      evaluate the slot \code{Map} of \code{RandVar} at \code{x}. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandVariable", x = "matrix", distr = "missing")}:
      evaluate the slot \code{Map} of \code{RandVar} at rows of \code{x}. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandVariable", x = "numeric", distr = "Distribution")}:
      evaluate the slot \code{Map} of \code{RandVar} at \code{x} assuming 
      a probability space with distribution \code{distr}. In case \code{x}
      does not lie in the support of \code{distr} \code{NA} is returned. }

    \item{evalRandVar}{\code{signature(RandVar = "EuclRandVariable", x = "matrix", distr = "Distribution")}:
      evaluate the slot \code{Map} of \code{RandVar} at rows of \code{x}
      assuming a probability space with distribution \code{distr}. For those
      rows of \code{x} which do not lie in the support of \code{distr} 
      \code{NA} is returned. }

    \item{imageDistr}{\code{signature(RandVar = "EuclRandVariable", distr = "Distribution")}:
      image distribution of \code{distr} under \code{RandVar}. Returns
      an object of class \code{"DistrList"}. }

    \item{dimension}{\code{signature(object = "EuclRandVariable")}:
      dimension of the Euclidean random variable. }

    \item{t}{\code{signature(x = "EuclRandVariable")}:
      returns an object of class \code{"EuclRandMatrix"} where the
      rhe results of the functions in the slot \code{Map} of \code{x} 
      are transposed. }

    \item{\%*\%}{\code{signature(x = "matrix", y = "EuclRandVariable")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}.   }

    \item{\%*\%}{\code{signature(x = "EuclRandVariable", y = "matrix")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{\%*\%}{\code{signature(x = "numeric", y = "EuclRandVariable")}: 
      generates an object of class \code{"EuclRandMatrix"} (1 x 1 matrix) 
      by multiplying (scalar/innner product) \code{x} and \code{y}. }

    \item{\%*\%}{\code{signature(x = "EuclRandVariable", y = "numeric")}: 
      generates an object of class \code{"EuclRandMatrix"} (1 x 1 matrix) 
      by multiplying (scalar/innner product) \code{x} and \code{y}. }

    \item{\%*\%}{\code{signature(x = "EuclRandVariable", y = "EuclRandVariable")}: 
      generates an object of class \code{"EuclRandMatrix"} (1 x 1 matrix) 
      by multiplying (scalar/innner product) \code{x} and \code{y}. }

    \item{\%*\%}{\code{signature(x = "EuclRandVariable", y = "EuclRandMatrix")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{\%*\%}{\code{signature(x = "EuclRandMatrix", y = "EuclRandVariable")}: 
      matrix multiplication of \code{x} and \code{y}. Generates 
      an object of class \code{"EuclRandMatrix"}. }

    \item{Arith}{\code{signature(e1 = "numeric", e2 = "EuclRandVariable")}:
      Given a numeric vector \code{e1}, a Euclidean random variable \code{e2} 
      and an arithmetic operator \code{op}, the Euclidean random variable 
      \code{e1 op e2} is returned. }

    \item{Arith}{\code{signature(e1 = "EuclRandVariable", e2 = "numeric")}:
      Given a numeric vector \code{e2}, a Euclidean random variable \code{e1} 
      and an arithmetic operator \code{op}, the Euclidean random variable 
      \code{e1 op e2} is returned. }

    \item{Arith}{\code{signature(e1 = "EuclRandVariable", e2 = "EuclRandVariable")}:
      Given two Euclidean random variables \code{e1}, \code{e2} and an 
      arithmetic operator \code{op}, the Euclidean random variable 
      \code{e1 op e2} is returned. }

    \item{Math}{\code{signature(x = "EuclRandVariable")}:
      Given a \code{"Math"} group generic \code{fct}, the Euclidean random 
      variable \code{fct(x)} is returned. }

    \item{E}{\code{signature(object = "UnivariateDistribution", fun = "EuclRandVariable", cond = "missing")}:
      expectation of \code{fun} under univariate distributions. }

    \item{E}{\code{signature(object = "AbscontDistribution", fun = "EuclRandVariable", cond = "missing")}:
      expectation of \code{fun} under absolutely continuous univariate distributions. }

    \item{E}{\code{signature(object = "DiscreteDistribution", fun = "EuclRandVariable", cond = "missing")}:
      expectation of \code{fun} under discrete univariate distributions. }

    \item{E}{\code{signature(object = "MultivariateDistribution", fun = "EuclRandVariable", cond = "missing")}: 
      expectation of \code{fun} under multivariate distributions. }

    \item{E}{\code{signature(object = "DiscreteMVDistribution", fun = "EuclRandVariable", cond = "missing")}: 
      expectation of \code{fun} under discrete multivariate distributions. }

    \item{E}{\code{signature(object = "UnivariateCondDistribution", fun = "EuclRandVariable", cond = "numeric")}: 
      conditional expectation of \code{fun} under conditional univariate distributions. }

    \item{E}{\code{signature(object = "UnivariateCondDistribution", fun = "EuclRandVariable", cond = "numeric")}: 
      conditional expectation of \code{fun} under absolutely continuous conditional univariate distributions. }

    \item{E}{\code{signature(object = "UnivariateCondDistribution", fun = "EuclRandVariable", cond = "numeric")}: 
      conditional expectation of \code{fun} under discrete conditional univariate distributions. }
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{EuclRandVariable}}, \code{\link{RandVariable-class}}, 
    \code{\link{EuclRandMatrix-class}}, \code{\link{EuclRandVarList-class}}, 
    \code{\link[distr]{Distribution-class}}, \code{\link[base:groupGeneric]{Arith}}, 
    \code{\link[base:groupGeneric]{Math}}, \code{\link[distrEx]{E}}}
\examples{
L1 <- list(function(x){x}, function(x){x^2}, function(x){x^3}, function(x){x^4})
L2 <- list(function(x){exp(x)}, function(x){abs(x)}, 
           function(x){sin(x)}, function(x){floor(x)})

R1 <- new("EuclRandVariable", Map = L1, Domain = Reals(), Range = Reals())
dimension(R1)
Map(R1)
Range(R1)
R1[2]
Map(R1[3])
Map(R1[c(1,2,4)])
Map(R1[2:4])
evalRandVar(R1, rnorm(1))
x <- as.matrix(rnorm(10))
res.R1 <- evalRandVar(R1, x)
res.R1[2,,] # results for Map(R1)[[2]](x)
res.R1[2,1,] # results for Map(R1)[[2]](x[1,])

R2 <- EuclRandVariable(L2, Domain = Reals(), dimension = 1)
dimension(R2)
DL1 <- imageDistr(R2, Norm())
plot(DL1)

Domain(R2) <- EuclideanSpace(dimension = 2)
Range(R2) <- EuclideanSpace(dimension = 2)
dimension(R2)
(X <- matrix(c(x, rnorm(10)), ncol = 2))
res.R2 <- evalRandVar(R2, X)
res.R2[3,,1] # results for Map(R2)[[3]](X[,1])

Map(log(abs(R2))) # "Math" group generic

# "Arith" group generic
Map(3 + R1)
Map(c(1,3,5) * R1)
try(1:5 * R1) # error
Map(1:2 * R2)
Map(R2 - 5)
Map(R1 ^ R1)

}
\concept{random variable}
\keyword{classes}
\keyword{arith}
\keyword{math}
