\name{get1Dderiv}
\alias{get1Dderiv}
\title{get1Dderiv}
\usage{
  get1Dderiv(w,A,B)
}

\description{
  The derivative of the objective function for the 1D-algorithm.
}

\arguments{
  \item{w}{A vector of length of \eqn{p}.}

  \item{A}{A \eqn{\sqrt{n}} estimate of an estimator's asymptotic 
    covariance matrix.}

  \item{B}{A \eqn{\sqrt{n}} estimate of the parameter associated with 
    the space we are enveloping. For our purposes this quantity 
    is either the outer product of the MLE of the mean-value 
    submodel parameter vector with itself or the outer product of the 
    MLE of the canonical submodel parameter vector with itself.}    
}

\value{
  \item{dF}{The value of the derivative of the objective function 
    for the 1D-algorithm evaluated at \code{w}, \code{A}, and \code{B}.}
}

\details{
  This function evaluates the derivative of the objective function 
  for the 1D-algorithm at \code{w}, \code{A}, and \code{B}. This is 
  needed in order to reliably find the maximum of the 1D-algorithm 
  objective function.
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
get1Dderiv(w = beta, A = avar, B = U)}
}
