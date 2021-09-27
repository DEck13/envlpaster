\name{get1Dobj}
\alias{get1Dobj}
\title{get1Dobj}
\usage{
  get1Dobj(w,A,B)
}

\description{
  The objective function for the 1D-algorithm.
}

\arguments{
  \item{w}{A vector of length of \code{p}.}

  \item{A}{A \eqn{\sqrt{n}} estimate of an estimator's asymptotic 
    covariance matrix.}

  \item{B}{A \eqn{\sqrt{n}} estimate of the parameter associated with 
    the space we are enveloping. for our purposes this quantity 
    is either the outer product of the MLE of the mean-value 
    submodel parameter vector with itself or the outer product of the 
    MLE of the canonical submodel parameter vector with itself.}    
}

\value{
  \item{Fw}{The value of the objective function for the 
  1D-algorithm evaluated at \code{w}, \code{A}, and \code{B}.}
}

\details{
  This function evaluates the objective function for the 
  1D-algorithm at \code{w}, \code{A}, and \code{B}. The maximizer 
  of this objective function is desired for a problem specific 
  \code{A} and \code{B}.
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
get1Dobj(w = beta, A = avar, B = U)}
}
