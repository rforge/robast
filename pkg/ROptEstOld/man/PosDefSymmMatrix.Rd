\name{PosDefSymmMatrix}
\alias{PosDefSymmMatrix}

\title{Generating function for PosDefSymmMatrix-class}
\description{Generates an object of class \code{"PosDefSymmMatrix"}.}
\usage{PosDefSymmMatrix(mat)}
\arguments{
  \item{mat}{A numeric positive-definite, symmetric 
    matrix with finite entries.}
}
\details{
  If \code{mat} is no matrix, \code{as.matrix} is applied.
}
\value{Object of class \code{"PosDefSymmMatrix"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{PosDefSymmMatrix-class}}}
\examples{
PosDefSymmMatrix(1)
PosDefSymmMatrix(diag(2))

## The function is currently defined as
function(mat){ 
    if(!is.matrix(mat)) mat <- as.matrix(mat)
    new("PosDefSymmMatrix", mat)
}
}
\concept{positive-definite matrix}
\concept{symmetric matrix}
\keyword{array}
