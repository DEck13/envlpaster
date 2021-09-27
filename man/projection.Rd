\name{projection}
\alias{projection}
\title{projection}
\usage{
  projection(a)
}

\arguments{
  \item{a}{An \eqn{n \times m} matrix}
}

\value{
  \item{pa}{The projection into the column space of \code{a}}
}

\description{
  Computes a projection into the column space of a 
  user specified matrix.
}

\examples{
\dontrun{library(envlpaster)
data(simdata30nodes)
data <- simdata30nodes.asterdata
nnode <- length(vars)
xnew <- as.matrix(simdata30nodes[,c(1:nnode)])
m1 <- aster(xnew, root, pred, fam, modmat)
avar <- m1$fisher
projection(a = avar)}
}
