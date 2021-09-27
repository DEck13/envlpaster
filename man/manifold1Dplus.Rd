\name{manifold1Dplus}
\alias{manifold1Dplus}
\title{manifold1Dplus}
\usage{
  manifold1Dplus(M,U,u)
}

\description{
  The 1D algorithm
}

\arguments{
  \item{M}{A \eqn{\sqrt{n}} estimate of an estimator's asymptotic 
    covariance matrix.}

  \item{U}{A \eqn{\sqrt{n}} estimate of the parameter associated with 
    the space we are enveloping. For our purposes this quantity 
    is either the outer product of the MLE of the mean-value 
    submodel parameter vector with itself or the outer product of the 
    MLE of the canonical submodel parameter vector with itself.}    

  \item{u}{The dimension of the envelope space assumed.}    
}

\value{
  \item{G}{A \eqn{\sqrt{n}} estimator of the basis matrix for the 
    envelope subspace. This matrix has \code{u} columns}
}

\details{
  This function calls \code{get1Dobj}, \code{get1Dini}, and \code{get1Dderiv} 
  in order to find
  \deqn{ \max_{w} \left[ \log(w^TMw) + \log(w^T(M+U)w) - 2\log(w^Tw) \right] }

  using Polak-Ribiere conjugate gradient in \code{optim}. This 
  maximization is conducted a total of \code{u} times and at each iteration 
  a vector belonging to the envelope space is returned. The vector 
  returned at a specific iteration is orthogonal to the vectors 
  returned at previous iterations. When finished, a basis matrix 
  for the envelope space is returned.  

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
manifold1Dplus(M = avar, U = U, u = 1)}
}
