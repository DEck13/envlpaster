\name{fit.boot}
\alias{fit.boot}
\title{fit.boot}
\usage{
  fit.boot(model, nboot, index, vectors = NULL, u = NULL,
    data, amat, newdata, modmat.new = NULL, renewdata = NULL, 
    fit.name = NULL, method = c("eigen","1d"), quiet = FALSE, m = 100)
}

\description{
  A parametric bootstrap procedure evaluated at an envelope estimator 
  of the submodel mean-value parameter vector \eqn{\tau} that was 
  obtained using reducing subspaces.  
}

\arguments{
  \item{model}{An aster model object.}

  \item{nboot}{The number of bootstrap iterations desired.}

  \item{index}{The indices denoting which components of the 
    canonical parameter vector are parameters of interest.}

  \item{vectors}{The indices denoting which reducing subspace 
    of Fisher information is desired to construct 
    envelope estimators.}

  \item{u}{The envelope model dimension.}

  \item{data}{An asterdata object corresponding to the 
    original data.}

  \item{amat}{This object can either be an array or a matrix. 
    It specifies a linear combination of mean-value parameters 
    that correspond to expected Darwinian fitness. See the 
    \code{aster} function help page in the original \code{aster} 
    package for more details.}

  \item{newdata}{A dataframe corresponding to hypothetical 
    individuals in which expected Darwinian fitness is to 
    be estimated.}    

  \item{modmat.new}{A model matrix corresponding to hypothetical 
    individuals in which expected Darwinian fitness is to 
    be estimated.}

  \item{renewdata}{A dataframe in long format corresponding to 
    hypothetical individuals in which expected Darwinian fitness 
    is to be estimated.}  

  \item{fit.name}{An expression that appears in the name of the 
    nodes that correspond to Darwinian fitness. This is only
    necessary if \code{renewdata} is not provided.}

  \item{method}{The procedure used to obtain envelope estimators.}    

  \item{quiet}{A logical argument. If FALSE, the function 
    displays how much time it takes to run \code{m} iterations.}

  \item{m}{The length of the output interval.}
}

\value{
  \item{u}{The dimension of the envelope space assumed}

  \item{table}{A table of output. The first two columns display 
    the envelope estimator of expected Darwinian fitness and its 
    bootstrapped standard error. The next two columns display the 
    MLE of expected Darwinian fitness and its bootstrapped 
    standard error. The last column displays the ratio of the 
    standard errors using MLE to those using envelope estimation. 
    Ratios greater than 1 indicate efficiency gains obtained 
    using envelope estimation.}

  \item{S}{The bootstrap estimator of the variability of 
    the partial envelope estimator.}

  \item{S2}{The bootstrap estimator of the variability of 
    the MLE.}

  \item{env.boot.out}{The realizations from the bootstrap procedure 
    using envelope methodology.}
    
  \item{MLE.boot.out}{The realizations from the bootstrap procedure 
    using maximum likelihood estimation.}
}

\details{
  This function implements the parametric bootstrap procedure given by 
  the algorithm presented below with respect to the mean-value 
  parameterization. This parametric bootstrap generates resamples from the distribution 
  evaluated at an envelope estimator of \eqn{\tau}. 

  The user specifies which vectors are used in order to construct envelope 
  estimators using the reducing subspace approach. The user also specifies which 
  method is to be used in order to calculate envelope estimators. When one 
  is using a partial envelope, then this function constructs envelope 
  estimators of \eqn{\upsilon} where we write \eqn{\tau} = \eqn{(\gamma^T,\upsilon^T)^T} 
  and \eqn{\upsilon} corresponds to aster model parameters of interest. 
  In applications, candidate reducing subspaces are indices of eigenvectors of \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}} 
  where \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}} is the part of \eqn{\hat{\Sigma}} 
  corresponding to our parameters of interest. These indices are specified 
  by \code{vectors}. When all of the components of \eqn{\tau} are components 
  of interest, then we write \eqn{\widehat{\Sigma}_{\upsilon,\upsilon} = \widehat{\Sigma}}. When data 
  is generated via the parametric bootstrap, it is the indices (not the 
  original reducing subspaces) that are used to construct envelope estimators 
  constructed using the generated data. The algorithm using reducing subspaces 
  is as follows: 

  \enumerate{
  \item[1.] Fit the aster model to the data and obtain 
    \eqn{\hat{\tau} = (\hat{\gamma}^T, \hat{\upsilon}^T)} and \eqn{\hat{\Sigma}} 
    from the aster model fit.
  \item[2.] Compute the envelope estimator of \eqn{\upsilon} in the original 
    sample, given as \eqn{\hat{\upsilon}_{env} = P_{\hat{G}}\hat{\upsilon}} 
    where \eqn{P_{\hat{G}}} is the projection into the reducing subspace of 
    \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}} specified by \code{vectors}. 
  \item[3.] Perform a parametric bootstrap by generating resamples from the 
  distribution evaluated at \eqn{\hat{\tau}_{env} = (\hat{\gamma}^T,\hat{\upsilon}_{env}^T)^T}. 
  For iteration \eqn{b=1,...,B} of the procedure: 
    \enumerate{
      \item[(3a)] Compute \eqn{\hat{\upsilon}^{(b)}} and 
        \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{(b)}} from the aster model fit to 
        the resampled data.
      \item[(3b)] Obtain \eqn{P_{\hat{G}}^{(b)}} as done in Step 2. 
      \item[(3c)] Compute \eqn{\hat{\upsilon}_{env}^{(b)} = P_{\hat{G}}^{(b)}\hat{\upsilon}^{(b)}} and 
        \eqn{\hat{\tau}_{env}^{(b)} = (\hat{\gamma}^{(b)^T},\hat{\upsilon}_{env}^{(b)^T})^T}.
      \item[(3d)] Store \eqn{h\left(\hat{\tau}_{env}^{(b)}\right)} where \eqn{h} 
        maps \eqn{\tau} to the parameterization of Darwinian fitness as determined 
        by \code{amat}.
    }
  }

  The parametric bootstrap procedure which uses the 1-d algorithm to construct 
  envelope estimators is analogous to the above algorithm. To use the 1-d 
  algorithm, the user specifies a candidate envelope model dimension \code{u} 
  and specifies \code{method = "1d"}. A parametric bootstrap generating resamples 
  from the distribution evaluated at the aster model MLE is also conducted by 
  this function.  
}

\references{
  Cook, R.D. and Zhang, X. (2014).
  Foundations for Envelope Models and Methods.
  \emph{JASA}, In Press.

  Cook, R.D. and Zhang, X. (2015).
  Algorithms for Envelope Estimation.
  \emph{Journal of Computational and Graphical Statistics}, 
  Published online. \doi{10.1080/10618600.2015.1029577}.

  Eck, D. J., Geyer, C. J., and Cook, R. D. (2016).
  Enveloping the aster model.
  \emph{in prep}.
}

\examples{
\dontrun{set.seed(13)
library(envlpaster)
library(aster2)
data(simdata30nodes)
data <- simdata30nodes.asterdata
nnode <- length(vars)
xnew <- as.matrix(simdata30nodes[,c(1:nnode)])
m1 <- aster(xnew, root, pred, fam, modmat)
target <- 5:9
indices <- c(1,2,4,5)
u <- length(indices)
nboot <- 2000; timer <- nboot/2
bar <- eigenboot(m1, nboot = nboot, index = target, 
  u = u, vectors = indices, data = data, m = timer)
bar}
}
