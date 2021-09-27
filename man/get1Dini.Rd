\name{get1Dini}
\alias{get1Dini}
\title{get1Dini}
\usage{
  get1Dini(A,B)
}

\description{
  Finds appropriate starting values for the 1D-algorithm.  
}

\arguments{
  \item{A}{A \eqn{\sqrt{n}} estimate of an estimator's asymptotic 
    covariance matrix.}

  \item{B}{A \eqn{\sqrt{n}} estimate of the parameter associated with 
    the space we are enveloping. For our purposes this quantity 
    is either the outer product of the MLE of the mean-value 
    submodel parameter vector with itself or the outer product of the 
    MLE of the canonical submodel parameter vector with itself.}    
}

\value{
  \item{w}{An appropriate starting value for the 1D-algorithm.}
}

\details{
  The 1D-algorithm is sensitive to starting values. Maximizers of 
  the 1D-algorithm objective function are close to eigenvectors of 
  \code{A} or (\code{A}+\code{B}). The vector, \code{w}, is the 
  eigenvector of \code{A} or (\code{A}+\code{B}) which produces the 
  largest value of the 1D-algorithm objective function. 
}

\references{
  Cook, R.D. and Zhang, X. (2014).
  Foundations for Envelope Models and Methods.
  \emph{JASA}, In Press.

  Cook, R.D. and Zhang, X. (2015).
  Algorithms for Envelope Estimation.
  \emph{Journal of Computational and Graphical Statistics}, 
  Published online. \doi{10.1080/10618600.2015.1029577}.  
}

\examples{
\dontrun{library(envlpaster)
data(simdata30nodes)
data <- simdata30nodes.asterdata
nnode <- length(vars)
xnew <- as.matrix(simdata30nodes[,c(1:nnode)])
m1 <- aster(xnew, root, pred, fam, modmat)
avar <- m1$fisher
beta <- m1$coef
U <- beta \%o\% beta
get1Dini(A = avar, B = U)}
}
