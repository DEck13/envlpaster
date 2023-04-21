\name{simdata30nodes}
\docType{data}
\alias{simdata30nodes}
\alias{variables}
\alias{root}
\alias{modmat}
\alias{formula}
\alias{xlevels}
\alias{terms}
\alias{simdata30nodes.asterdata}
\title{A generated aster data set with 30 nodes}
\description{
  Simulated data for an aster analysis. Loads 7 objects.
}

\usage{data(simdata30nodes)}
\format{The data frame with records for 3000 organisms over 10 years. 
  The dataset corresponding to our aster analysis. The following four descriptions explain the elements of this dataset.

  \describe{

    \item{u}{Indicates survival for each of the 10 years.}
    \item{w}{Counts offspring for each of the 10 years.}
    \item{v}{Indicates if \code{w > 0} for each of the 10 years.}
    \item{z}{A covariate of potential interest, 10 in total.}
    \item{variables}{Character vector giving the names of the variables in 
      the graph.}
    \item{root}{The root data. For \code{aster.default} an 
      \code{nind} by \code{nnode} matrix, for \code{aster.formula} 
      an \code{nind * nnode} vector.}
    \item{modmat}{An \code{nind} by \code{nnode} by \code{ncoef} 
      three-dimensional array, the model matrix. \code{aster.formula} 
      constructs such a modmat from its formula, the data frame data,
      and the variables in the environment of the formula.}
    \item{formula}{Necessary for changing to class \code{aster.formula}.}
    \item{xlevels}{Necessary for changing to class \code{aster.formula}.}
    \item{terms}{Necessary for changing to class \code{aster.formula}.}
    \item{simdata30nodes.asterdata}{An object of class \code{asterdata} 
      corresponding to \code{simdata30nodes}.}
  }
}

\details{
  This object contains an aster data set in wide form, an 
  object of class \code{asterdata} corresponding to the 
  original data set, and vectors specifying the graphical 
  structure of the aster model.

  There are 3000 simulated individuals in this aster analysis. 
  Our data is generated in two parts. The first part follows 
  Technical report 671 (TR 671) on Charlie Geyer's Aster 
  Models for Life History Analysis	webpage. For our data, 
  \code{nind = 3000}, \code{ntime = 10}, \code{psurv = 0.95}, \code{prepr = 0.7}, \code{mpois = 1}, 
  and the seed is set at \code{set.seed(13)} which is different 
  from the original simulation setup. 

  We follow the model construction in TR 671 through \code{out6}. 
  We then generate a new dataset from the aster model where 
  the components of the submodel mean-value parameter vector 
  \eqn{\tau} corresponding to Darwinian fitness is in the 
  space spanned by the first, second, and fourth eigenvectors 
  of Fisher information.
}


%\source{}
\references{
  Geyer, C. J. and Shaw, R. G. (2009). 
  Model Selection in Estimation of Fitness Landscapes. Technical Report No. 671. School of Statistics, University of Minnesota. 
  \url{http://conservancy.umn.edu/handle/11299/56219}.
}

\keyword{datasets}