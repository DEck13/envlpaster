\name{scanner}
\alias{scanner}
\title{scanner}
\usage{
  scanner(M, coef, u)
}

\description{
  A diagnostic assessing the potential benefits of envelope 
  methodology for a particular aster model.
}

\arguments{
  \item{M}{A \eqn{\sqrt{n}} estimate of an estimator's asymptotic 
    covariance matrix.}

  \item{coef}{The MLE of the parameter of interest.}

  \item{u}{The dimension of the envelope space assumed.}
}

\value{
  \item{indices}{The indices of the \code{u} most relevant 
    eigenvectors of \code{M} to the construction of \code{coef}.}

  \item{table}{An output table. The first column is the projection 
    of \code{coef} into the space spanned by the eigenvectors given 
    by the \code{indices}.}

  \item{G}{The \code{u} most relevant eigenvectors of \code{M} 
    to the construction of \code{coef}.}

  \item{prop}{the sum of the remaining eigenvalues of \code{M} divided 
    by the sum of all of the eigenvalues of \code{M}.}
}

\details{
  This function provides users with a rough diagnostic for the performance 
  of an envelope estimator at a specific dimension. We can see how close 
  a particular potential envelope estimator is to the MLE as well as the 
  proportion of variation that will be discarded when using envelope 
  estimation. This amount of variation discarded is optimistic since it 
  does not account for variability associated with estimating the 
  projection into the envelope space.
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
scanner(M = avar, coef = beta, u = 1)}
}
