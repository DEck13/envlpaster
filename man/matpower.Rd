\name{matpower}
\alias{matpower}
\title{matpower}
\usage{
  matpower(a,alpha)
}

\description{
  Raises a matrix to a power.
}

\arguments{
  \item{a}{An \eqn{n \times n} matrix}

  \item{alpha}{The power which to raise the matrix to.}
}

\value{
  \item{ai}{The matrix raised to the desired power}
}


\examples{
\dontrun{library(envlpaster)
data(simdata30nodes)
data <- simdata30nodes.asterdata
nnode <- length(vars)
xnew <- as.matrix(simdata30nodes[,c(1:nnode)])
m1 <- aster(xnew, root, pred, fam, modmat)
avar <- m1$fisher
matpower(a = avar, alpha = -1)}
}